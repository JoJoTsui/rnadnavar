from pathlib import Path

from pptx import Presentation
from pptx.util import Inches, Pt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs/presentation/FASTQ_TO_RESCUED_VCF_OVERVIEW.pptx"
IMG = ROOT / "docs/images/rnadnavar_schemav3.png"

prs = Presentation()


def add_title_slide(title: str, subtitle: str) -> None:
    slide = prs.slides.add_slide(prs.slide_layouts[0])
    slide.shapes.title.text = title
    slide.placeholders[1].text = subtitle


def add_bullets_slide(title: str, bullets: list[str]) -> None:
    slide = prs.slides.add_slide(prs.slide_layouts[1])
    slide.shapes.title.text = title
    tf = slide.shapes.placeholders[1].text_frame
    tf.clear()
    for i, line in enumerate(bullets):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        p.text = line
        p.level = 0
        p.font.size = Pt(22)


def add_workflow_image_slide() -> None:
    slide = prs.slides.add_slide(prs.slide_layouts[5])
    slide.shapes.title.text = "Pipeline Workflow at a Glance"
    if IMG.exists():
        slide.shapes.add_picture(str(IMG), Inches(0.6), Inches(1.2), width=Inches(12.0))
    textbox = slide.shapes.add_textbox(
        Inches(0.7), Inches(6.2), Inches(12.0), Inches(0.8)
    )
    tf = textbox.text_frame
    tf.text = "DNA and RNA paths run in parallel, then cross over in rescue; optional RNA realignment enables second rescue."
    tf.paragraphs[0].font.size = Pt(18)


add_title_slide(
    "rnadnavar: FASTQ to Rescued VCF",
    "Code-accurate documentation update | Conference overview",
)

add_bullets_slide(
    "Why pair DNA and RNA",
    [
        "DNA and RNA provide complementary evidence for somatic variants",
        "Consensus reduces single-caller noise",
        "Rescue keeps cross-modality supported variants",
    ],
)

add_bullets_slide(
    "Input Labeling in This Run",
    [
        "status 0: COO8801DN (DNA normal)",
        "status 1: COO8801DT (DNA tumor)",
        "status 2: COO8801RT (RNA tumor)",
        "patient id COO8801 links all rows in one case set",
    ],
)

add_workflow_image_slide()

add_bullets_slide(
    "Consensus Logic (within modality)",
    [
        "Aggregate caller VCFs by normalized variant key",
        "Apply thresholds: SNV 2, indel 2 (this run)",
        "Clear majority -> unified class",
        "Tie at top class -> Artifact",
    ],
)

add_bullets_slide(
    "Rescue Logic (cross modality)",
    [
        "Combine DNA and RNA consensus with individual caller support",
        "Agreement keeps class; disagreement handled with support-aware rules",
        "First rescue uses original RNA consensus",
        "Second rescue can use realigned RNA consensus",
    ],
)

add_bullets_slide(
    "FILTER View for Presentation",
    [
        "Somatic: prioritized candidates",
        "Germline: likely inherited",
        "Reference: evidence favors non-variant/reference state",
        "Detailed classes Artifact/NoConsensus/RNAedit are in the full guide",
    ],
)

add_bullets_slide(
    "Verified Outputs and Run Settings",
    [
        "Consensus and rescue files verified in output/COO8801.shared",
        "realignment_mode = vcf",
        "tools include consensus,rescue,realignment,vep",
        "No algorithm or code changes in this update",
    ],
)

prs.save(str(OUT))
print(f"Wrote {OUT}")
