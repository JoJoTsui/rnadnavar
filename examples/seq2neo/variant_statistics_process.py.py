import os
import glob
import gzip
from cyvcf2 import VCF
import pandas as pd
from concurrent.futures import ProcessPoolExecutor


def find_vcf_files(base_dir, pattern):
    """
    Recursively finds VCF files matching a given pattern within a base directory.
    """
    search_path = os.path.join(base_dir, '**', pattern)
    return glob.glob(search_path, recursive=True)


def parse_main_vcf(vcf_file_path):
    """
    Parses a main VCF file and extracts CHROM, POS, REF, ALT, FILTER,
    CALLERS_SUPPORT, COSMIC_ID, and GNOMAD_AF information using cyvcf2.
    """
    records_data = []
    try:
        vcf_reader = VCF(vcf_file_path)
        for record in vcf_reader:
            callers_support = record.INFO.get('CALLERS_SUPPORT')
            cosmic_id = record.INFO.get('COSMIC_ID')
            gnomad_af = record.INFO.get('GNOMAD_AF')

            records_data.append({
                'CHROM': record.CHROM,
                'POS': record.POS,
                'REF': record.REF,
                'ALT': record.ALT[0] if record.ALT else None,
                'FILTER': record.FILTERS[0] if record.FILTERS else 'PASS',
                'CALLERS_SUPPORT': callers_support,
                'COSMIC_ID': cosmic_id if cosmic_id is not None else None,
                'GNOMAD_AF': gnomad_af if gnomad_af is not None else None
            })
    except Exception as e:
        print(f"Error parsing {vcf_file_path}: {e}")

    records_df = pd.DataFrame(records_data)
    print(f"records_df size", records_df.shape)
    return records_df


def process_single_caller(chroms, poss, mask, vcf_path, sample_suffix):
    """
    处理单个caller的VCF：顺序扫描构建字典，然后按位置查表。
    返回两个等长列表: dp_list, ad_list
    """
    n = len(chroms)
    dp_list = [None] * n
    ad_list = [None] * n

    try:
        vcf_reader = VCF(vcf_path)
        # 查找匹配后缀的样本索引
        sample_index = -1
        for i, name in enumerate(vcf_reader.samples):
            if name.endswith(sample_suffix):
                sample_index = i
                break
        if sample_index == -1:
            print(f"Warning: No sample with suffix '{sample_suffix}' in {vcf_path}")
            vcf_reader.close()
            return dp_list, ad_list

        # 1) 顺序扫描，构建位置 → (DP, AD) 的字典
        pos_dict = {}
        for record in vcf_reader:
            dp = record.format('DP')[sample_index][0] if 'DP' in record.FORMAT else None
            ad = record.format('AD')[sample_index] if 'AD' in record.FORMAT else None
            ad_str = ','.join(str(a) for a in ad) if ad is not None else None
            pos_dict[(record.CHROM, record.POS)] = (dp, ad_str)
        vcf_reader.close()

        # 2) 根据主VCF的位置及mask批量查表
        for i in range(n):
            if mask[i]:
                key = (chroms[i], poss[i])
                if key in pos_dict:
                    dp_list[i], ad_list[i] = pos_dict[key]

    except Exception as e:
        print(f"Error processing {vcf_path}: {e}")

    return dp_list, ad_list


def process_callers_support(df, base_output_dir, project_name):
    """
    并行读取6个caller VCF，将DP/AD信息填入DataFrame。
    """
    # 定义各caller的路径信息
    prefix = project_name
    callers_info = {
        'DNA_deepsomatic': {
            'path': os.path.join(base_output_dir, project_name,
                                 "variant_calling/deepsomatic",
                                 f"{prefix}DT_vs_{prefix}DN"),
            'pattern': "*.deepsomatic.vcf.gz",
            'sample_suffix': 'DT'
        },
        'DNA_mutect2': {
            'path': os.path.join(base_output_dir, project_name,
                                 "variant_calling/mutect2",
                                 f"{prefix}DT_vs_{prefix}DN"),
            'pattern': "*.mutect2.filtered.vcf.gz",
            'sample_suffix': 'DT'
        },
        'DNA_strelka': {
            'path': os.path.join(base_output_dir, project_name,
                                 "variant_calling/strelka",
                                 f"{prefix}DT_vs_{prefix}DN"),
            'pattern': "*.strelka.variants.vcf.gz",
            'sample_suffix': 'TUMOR'
        },
        'RNA_deepsomatic': {
            'path': os.path.join(base_output_dir, project_name,
                                 "vcf_realignment/variant_calling/deepsomatic",
                                 f"{prefix}RT_realign_vs_{prefix}DN"),
            'pattern': "*.deepsomatic.vcf.gz",
            'sample_suffix': 'RT'
        },
        'RNA_mutect2': {
            'path': os.path.join(base_output_dir, project_name,
                                 "vcf_realignment/variant_calling/mutect2",
                                 f"{prefix}RT_realign_vs_{prefix}DN"),
            'pattern': "*.mutect2.filtered.vcf.gz",
            'sample_suffix': 'RT'
        },
        'RNA_strelka': {
            'path': os.path.join(base_output_dir, project_name,
                                 "vcf_realignment/variant_calling/strelka",
                                 f"{prefix}RT_realign_vs_{prefix}DN"),
            'pattern': "*.strelka.variants.vcf.gz",
            'sample_suffix': 'TUMOR'
        }
    }

    # 1) 预处理：提取所有任务信息（找到VCF文件、生成mask）
    tasks = []
    for caller, info in callers_info.items():
        vcf_files = find_vcf_files(info['path'], info['pattern'])
        if not vcf_files:
            print(f"No VCF files found for {caller} in {info['path']}")
            continue
        vcf_path = vcf_files[0]   # 假设每个caller只有一个VCF

        # 创建布尔mask：该位点的CALLERS_SUPPORT中是否包含当前caller
        mask = df['CALLERS_SUPPORT'].apply(
            lambda x, c=caller: c in x if isinstance(x, str) else False
        ).tolist()

        tasks.append((caller, vcf_path, info['sample_suffix'], mask))

    # 2) 提取所有caller共用的染色体和位置列表（仅传递一次）
    chroms = df['CHROM'].tolist()
    poss   = df['POS'].tolist()

    # 3) 初始化新列（避免后续赋值时键不存在）
    for caller, _, _, _ in tasks:
        df[f'{caller}_DP'] = None
        df[f'{caller}_AD'] = None

    # 4) 多进程并行处理各个caller
    with ProcessPoolExecutor(max_workers=min(6, len(tasks))) as executor:
        futures = {}
        for caller, vcf_path, sample_suffix, mask in tasks:
            future = executor.submit(
                process_single_caller,
                chroms, poss, mask, vcf_path, sample_suffix
            )
            futures[future] = caller

        # 收集结果并回填DataFrame
        for future in futures:
            caller = futures[future]
            try:
                dp_list, ad_list = future.result()
                df[f'{caller}_DP'] = dp_list
                df[f'{caller}_AD'] = ad_list
            except Exception as e:
                print(f"Error retrieving results for {caller}: {e}")

    return df


def save_dataframe_to_csv(df, output_csv_path):
    """保存DataFrame为CSV文件"""
    try:
        df.to_csv(output_csv_path, index=False)
        print(f"Data successfully saved to {output_csv_path}")
    except Exception as e:
        print(f"Error saving DataFrame to CSV at {output_csv_path}: {e}")


def process_project_vcf_data(base_output_dir, project_name):
    """处理单个项目的完整流程"""
    print(f"Processing project: {project_name}")

    # 1. 查找主VCF文件
    prefix = project_name
    main_vcf_dir = os.path.join(
        base_output_dir, project_name,
        "vcf_realignment/rescue",
        f"{prefix}DT_vs_{prefix}DN_rescued_{prefix}RT_realign_vs_{prefix}DN"
    )
    main_vcf_pattern = "*.filtered.vcf.stripped.vcf.gz"
    main_vcf_files = find_vcf_files(main_vcf_dir, main_vcf_pattern)

    if not main_vcf_files:
        print(f"No main VCF file found for {project_name} in {main_vcf_dir}")
        return

    main_vcf_file_path = main_vcf_files[0]
    print(f"Found main VCF file: {main_vcf_file_path}")

    # 2. 解析主VCF
    df = parse_main_vcf(main_vcf_file_path)
    if df.empty:
        print(f"No records parsed from main VCF file: {main_vcf_file_path}")
        return

    # 3. 并行提取所有caller的DP/AD信息
    df = process_callers_support(df, base_output_dir, project_name)

    # 4. 保存结果
    output_csv_path = os.path.join(base_output_dir, f"{project_name}_processed_vcf_data.csv")
    save_dataframe_to_csv(df, output_csv_path)
    print(f"Finished processing project: {project_name}")


if __name__ == '__main__':
    base_output_dir = r"/t9k/mnt/WorkSpace/data/ngs/liuxin/seq2neo/output"

    # 自动获取 base_output_dir 下所有子目录作为项目名
    projects_to_process = [
        d for d in os.listdir(base_output_dir)
        if os.path.isdir(os.path.join(base_output_dir, d)) and not d.startswith('.')
    ]

    failed_projects = []  # 新增：记录失败项目

    for project_name in projects_to_process:
        try:
            process_project_vcf_data(base_output_dir, project_name)
        except Exception as e:
            error_msg = f"Project '{project_name}' failed with error: {e}"
            print(error_msg)
            failed_projects.append((project_name, str(e)))

    # 输出失败项目汇总
    if failed_projects:
        print("\n=== Summary of Failed Projects ===")
        for proj, err in failed_projects:
            print(f"  - {proj}: {err}")
    else:
        print("\nAll projects processed successfully.")