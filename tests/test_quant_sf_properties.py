"""
Property-based tests for quant.sf Column Count.

# Feature: neoantigen-workflow, Property 10: quant.sf Column Count

Validates: Requirements 4a.4, 8.4

For any quant.sf file produced by SALMON_QUANT, every non-header line must
contain exactly 5 tab-separated columns: Name, Length, EffectiveLength, TPM,
NumReads.

The quant.sf format (Salmon v1.11.4 spec):
- Header line: Name\tLength\tEffectiveLength\tTPM\tNumReads
- Data lines:  <transcript_id>\t<length>\t<eff_length>\t<tpm>\t<num_reads>
"""

from hypothesis import given, settings
from hypothesis import strategies as st

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

QUANT_SF_HEADER = "Name\tLength\tEffectiveLength\tTPM\tNumReads"
QUANT_SF_COLUMNS = ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]

# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

# Transcript name: plain or Gencode pipe-delimited format
transcript_name_st = st.one_of(
    # Plain transcript ID (e.g. ENST00000456328.2)
    st.text(
        alphabet=st.characters(whitelist_categories=("Lu", "Ll", "Nd"), whitelist_characters="._-"),
        min_size=1,
        max_size=30,
    ),
    # Gencode pipe-delimited format (e.g. ENST00000456328.2|ENSG00000223972.5|...)
    st.builds(
        lambda parts: "|".join(parts),
        st.lists(
            st.text(
                alphabet=st.characters(whitelist_categories=("Lu", "Ll", "Nd"), whitelist_characters="._-"),
                min_size=1,
                max_size=20,
            ),
            min_size=2,
            max_size=5,
        ),
    ),
)

# A single quant.sf data line with exactly 5 tab-separated fields
data_line_st = st.builds(
    lambda name, length, eff_length, tpm, num_reads: (
        f"{name}\t{length}\t{eff_length:.6f}\t{tpm:.6f}\t{num_reads:.6f}"
    ),
    name=transcript_name_st,
    length=st.integers(min_value=100, max_value=10000),
    eff_length=st.floats(min_value=1.0, max_value=10000.0, allow_nan=False, allow_infinity=False),
    tpm=st.floats(min_value=0.0, max_value=1e6, allow_nan=False, allow_infinity=False),
    num_reads=st.floats(min_value=0.0, max_value=1e8, allow_nan=False, allow_infinity=False),
)

# A complete quant.sf file content (header + 1-100 data lines)
quant_sf_content_st = st.builds(
    lambda lines: QUANT_SF_HEADER + "\n" + "\n".join(lines),
    lines=st.lists(data_line_st, min_size=1, max_size=100),
)


# ---------------------------------------------------------------------------
# Helper: parse quant.sf content
# ---------------------------------------------------------------------------


def parse_quant_sf(content: str) -> tuple[str, list[str]]:
    """
    Parse quant.sf content into (header_line, data_lines).
    Strips trailing whitespace/newlines from each line.
    """
    lines = [line for line in content.splitlines() if line.strip()]
    header = lines[0] if lines else ""
    data_lines = lines[1:] if len(lines) > 1 else []
    return header, data_lines


# ---------------------------------------------------------------------------
# Property 10: quant.sf Column Count — data lines have exactly 5 columns
# ---------------------------------------------------------------------------


@given(content=quant_sf_content_st)
@settings(max_examples=100)
def test_p10_every_data_line_has_exactly_5_columns(content):
    """
    **Validates: Requirements 4a.4, 8.4**

    Property 10: quant.sf Column Count.
    Every non-header line in a quant.sf file must have exactly 5
    tab-separated columns.

    # Feature: neoantigen-workflow, Property 10: quant.sf Column Count
    """
    _, data_lines = parse_quant_sf(content)

    assert len(data_lines) >= 1, "quant.sf must have at least one data line"

    for i, line in enumerate(data_lines):
        columns = line.split("\t")
        assert len(columns) == 5, (
            f"Data line {i + 1} has {len(columns)} column(s), expected exactly 5. "
            f"Line content: {line!r}"
        )


@given(content=quant_sf_content_st)
@settings(max_examples=100)
def test_p10_header_line_has_exactly_5_columns_with_correct_names(content):
    """
    **Validates: Requirements 4a.4, 8.4**

    Property 10: quant.sf Column Count — header validation.
    The header line must have exactly 5 tab-separated columns with the
    correct canonical names: Name, Length, EffectiveLength, TPM, NumReads.

    # Feature: neoantigen-workflow, Property 10: quant.sf Column Count
    """
    header, _ = parse_quant_sf(content)

    columns = header.split("\t")
    assert len(columns) == 5, (
        f"Header has {len(columns)} column(s), expected exactly 5. "
        f"Header: {header!r}"
    )
    assert columns == QUANT_SF_COLUMNS, (
        f"Header columns {columns} do not match expected {QUANT_SF_COLUMNS}."
    )
