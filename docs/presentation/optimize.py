#!/usr/bin/env python3
"""
Generates an optimized, 4-slide version of the FASTQ to Rescued VCF Overview.
It includes a clean layout, a programmatic flowchart, and optimized Chinese presentation notes.

Requirements:
    pip install python-pptx
"""

import os
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE

def set_elegant_title(title_shape, text):
    """Applies a consistent, elegant styling to title shapes."""
    title_shape.text = text
    for paragraph in title_shape.text_frame.paragraphs:
        paragraph.alignment = PP_ALIGN.LEFT
        paragraph.font.name = 'Arial'
        paragraph.font.size = Pt(36)
        paragraph.font.color.rgb = RGBColor(0x00, 0x33, 0x66) # Dark Blue

def add_bullet_points(body_shape, points):
    """Adds cleanly formatted bullet points."""
    text_frame = body_shape.text_frame
    text_frame.clear()
    for i, (point_text, level) in enumerate(points):
        p = text_frame.add_paragraph() if i > 0 else text_frame.paragraphs[0]
        p.text = point_text
        p.level = level
        p.font.name = 'Arial'
        p.font.size = Pt(22 - (level * 4)) # Scale down per level
        p.font.color.rgb = RGBColor(0x33, 0x33, 0x33) # Dark Gray

def add_flowchart(slide):
    """Draws a vertical workflow flowchart using PPTX shapes."""
    steps = [
        "1. Input FASTQ\n(Status 0, 1, 2)",
        "2. Alignment & Preprocessing\n(DNA/RNA Parallel Tracks)",
        "3. Variant Calling\n(DeepSomatic, Mutect2, Strelka)",
        "4. Within-Modality Consensus\n(Voting & Filtering)",
        "5. Cross-Modality Rescue\n(DNA & RNA Integration)",
        "6. Final VCF Labeled FILTER\n(Somatic, Germline, Reference)"
    ]
    
    start_top = Inches(1.5)
    left_margin = Inches(2.5)
    box_width = Inches(5.0)
    box_height = Inches(0.8)
    vertical_spacing = Inches(1.1)

    for i, step_text in enumerate(steps):
        # Draw Box
        top = start_top + (i * vertical_spacing)
        shape = slide.shapes.add_shape(
            MSO_SHAPE.ROUNDED_RECTANGLE, left_margin, top, box_width, box_height
        )
        
        # Style Box
        shape.fill.solid()
        shape.fill.fore_color.rgb = RGBColor(0xE6, 0xF0, 0xFA) # Light Blue
        shape.line.color.rgb = RGBColor(0x00, 0x4C, 0x99) # Border Blue
        shape.line.width = Pt(1.5)
        
        # Add Text
        text_frame = shape.text_frame
        text_frame.text = step_text
        for p in text_frame.paragraphs:
            p.alignment = PP_ALIGN.CENTER
            p.font.name = 'Arial'
            p.font.size = Pt(16)
            p.font.bold = True
            p.font.color.rgb = RGBColor(0x00, 0x33, 0x66)
            
        # Draw Arrow down to next box (except for last one)
        if i < len(steps) - 1:
            arrow_top = top + box_height
            arrow_height = vertical_spacing - box_height
            arrow = slide.shapes.add_shape(
                MSO_SHAPE.DOWN_ARROW, left_margin + box_width/2 - Inches(0.15), 
                arrow_top + Inches(0.05), Inches(0.3), arrow_height - Inches(0.1)
            )
            arrow.fill.solid()
            arrow.fill.fore_color.rgb = RGBColor(0x99, 0x99, 0x99)
            arrow.line.fill.background()

def create_optimized_presentation():
    # If the original template exists, load it to keep the master slides/theme
    # Otherwise, start fresh. We append slides to not break existing layout indices.
    template_path = "docs/presentation/FASTQ_TO_RESCUED_VCF_OVERVIEW.pptx"
    
    if os.path.exists(template_path):
        prs = Presentation(template_path)
    else:
        prs = Presentation()

    title_layout = prs.slide_layouts[0]
    content_layout = prs.slide_layouts[1]

    # --- SLIDE 1: Title Slide ---
    slide1 = prs.slides.add_slide(title_layout)
    title = slide1.shapes.title
    subtitle = slide1.placeholders[1]
    
    title.text = "rnadnavar: From FASTQ to Labeled FILTER"
    subtitle.text = "Optimized Variant Classification Workflow\nEnd-to-End Tracking and Rescue Integration"
    
    slide1.notes_slide.notes_text_frame.text = (
        "欢迎了解 rnadnavar 分析流程。\n"
        "本简报将带您鸟瞰从原始测序数据（FASTQ）到最终变异标记（Labeled VCF FILTER）的全过程。\n"
        "我们将把原始的详细步骤精简为核心的三大环节：变异检测原理、组学内部共识（Consensus）与跨组学拯救（Rescue），为您提供清晰的全景视图。"
    )

    # --- SLIDE 2: Workflow Overview (Flowchart) ---
    slide2 = prs.slides.add_slide(content_layout)
    set_elegant_title(slide2.shapes.title, "End-to-End Workflow Architecture")
    # Remove default text box to make room for flowchart
    if len(slide2.placeholders) > 1:
        sp = slide2.placeholders[1]
        sp.element.getparent().remove(sp.element)
    
    add_flowchart(slide2)
    
    slide2.notes_slide.notes_text_frame.text = (
        "如图所示，这是整个数据处理的完整流程图。\n"
        "1. 数据流分为并行处理的 DNA 和 RNA 两条主线。\n"
        "2. 数据经过高质量的比对（Alignment）后，交由多种检测器（如 DeepSomatic）进行变异检测（Variant Calling）。\n"
        "3. 接着，基于投票机制生成单组学共识（Modality Consensus）。\n"
        "4. 最后，在 Rescue 阶段进行跨组学相互印证和融合，输出最终的高可信 VCF 过滤标签。"
    )

    # --- SLIDE 3: Caller Classification & Consensus ---
    slide3 = prs.slides.add_slide(content_layout)
    set_elegant_title(slide3.shapes.title, "Variant Calling & Consensus Logic")
    
    points_s3 = [
        ("Per-Caller Classification Mapping:", 0),
        ("DeepSomatic & Mutect2 map direct filters to Somatic/Germline/Artifact.", 1),
        ("Strelka uses normal tier (NT) and read depths to derive biological class.", 1),
        ("Within-Modality Consensus (Voting):", 0),
        ("Variants are grouped by normalized keys separately in DNA/RNA tracks.", 1),
        ("Strict thresholds applied (e.g., SNV ≥ 2, Indel ≥ 2 valid caller support).", 1),
        ("A clear majority establishes the label; ties default to 'Artifact'.", 1)
    ]
    add_bullet_points(slide3.placeholders[1], points_s3)
    
    slide3.notes_slide.notes_text_frame.text = (
        "本页详细说明了单组学内部的数据处理逻辑。\n"
        "首先，各检测软件（如 DeepSomatic, Mutect2 和 Strelka）根据其独有的算法生成基础的突变分类标签。\n"
        "接着，系统会对这些标签进行汇总，即“共识（Consensus）”阶段：系统基于同类归一化后的变异进行多数投票。\n"
        "只有满足最低支持软件数量（例如至少2个软件支持）的突变才会获得有效共识标签。若各软件意见平局或不明确，为了严谨起见，通常会被保守标记为 Artifact（假阳性/伪影）。"
    )

    # --- SLIDE 4: Cross-Modality Rescue & Final Output ---
    slide4 = prs.slides.add_slide(content_layout)
    set_elegant_title(slide4.shapes.title, "Cross-Modality Rescue & Final Categories")
    
    points_s4 = [
        ("Cross-Modality Integration (Rescue):", 0),
        ("Combines DNA and RNA consensus labels on a per-variant basis.", 1),
        ("Agreeing labels validate the variant strongly.", 1),
        ("Conflicting labels undergo a support-aware resolution algorithm.", 1),
        ("Optional second rescue phase uses newly realigned RNA data.", 1),
        ("Final VCF FILTER Categories Delivered:", 0),
        ("Somatic: High-confidence acquired tumor mutations.", 1),
        ("Germline: Highly likely inherited variants.", 1),
        ("Reference: Insufficient variant evidence or reference-like state.", 1)
    ]
    add_bullet_points(slide4.placeholders[1], points_s4)
    
    slide4.notes_slide.notes_text_frame.text = (
        "跨组学拯救（Rescue）是流程的最高决断阶段。\n"
        "在此阶段，系统比对 DNA 与 RNA 的共识结果。如果双方一致，突变可信度极高；如果不一致，算法会基于多维度证据支持度进行高级决策（可结合可选的 RNA 二次比对进行补充拯救）。\n"
        "最终，变异结果会被清晰、精准地分类输出在 VCF 文件的 FILTER 列中：\n"
        "1. Somatic (体细胞突变)\n"
        "2. Germline (生殖细胞突变)\n"
        "3. Reference (参考集/非突变)\n"
        "极大降低了人工复核的复杂度。"
    )

    # Save the output
    output_path = "docs/presentation/OPTIMIZED_FASTQ_TO_RESCUED_VCF.pptx"
    prs.save(output_path)
    print(f"✅ Successfully generated optimized presentation: {output_path}")
    print("If you appended to the original file, you can now safely delete the older original slides inside the presentation.")

if __name__ == "__main__":
    create_optimized_presentation()