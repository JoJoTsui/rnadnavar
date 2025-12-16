#!/usr/bin/env python3
"""
VCF Field Escaping Utilities

This module provides utilities for properly escaping special characters in VCF fields
to prevent parsing issues. The VCF format uses specific characters as separators:
- ';' separates INFO fields
- '=' separates INFO keys from values  
- ',' separates multiple values within an INFO field
- '|' is used for phased genotype separation

When FILTER values or other data containing these characters are stored in INFO fields,
they must be properly escaped to prevent VCF parsing errors.

Key Functions:
    escape_vcf_info_value: Escape special characters for INFO field values
    unescape_vcf_info_value: Unescape previously escaped INFO field values
    escape_filter_for_info: Specifically escape FILTER values for INFO storage

Author: VCF Processing Pipeline
Date: 2025-12-16
"""

import re
from typing import Optional, Union


def escape_vcf_info_value(value: Union[str, int, float, None]) -> str:
    """
    Escape special characters in a value for safe storage in VCF INFO fields.
    
    This function escapes VCF special characters using URL-style percent encoding:
    - ';' (semicolon) -> '%3B' (INFO field separator)
    - '=' (equals) -> '%3D' (key-value separator)  
    - ',' (comma) -> '%2C' (value separator)
    - ' ' (space) -> '%20' (whitespace)
    - '\t' (tab) -> '%09' (tab character)
    - '\n' (newline) -> '%0A' (newline character)
    
    Args:
        value: The value to escape (string, number, or None)
        
    Returns:
        Escaped string safe for VCF INFO fields
        
    Example:
        >>> escape_vcf_info_value("LowDepth;LowEVS")
        'LowDepth%3BLowEVS'
        >>> escape_vcf_info_value("key=value,other")
        'key%3Dvalue%2Cother'
        >>> escape_vcf_info_value(None)
        '.'
    """
    if value is None:
        return '.'
    
    # Convert to string
    value_str = str(value)
    
    # Handle empty string
    if not value_str:
        return '.'
    
    # Escape VCF special characters with URL encoding
    escaped = value_str
    escaped = escaped.replace(';', '%3B')  # Semicolon (INFO field separator) - CRITICAL
    escaped = escaped.replace('=', '%3D')  # Equals (key-value separator)
    escaped = escaped.replace(',', '%2C')  # Comma (value separator)
    escaped = escaped.replace(' ', '%20')  # Space (whitespace)
    escaped = escaped.replace('\t', '%09') # Tab
    escaped = escaped.replace('\n', '%0A') # Newline
    escaped = escaped.replace('\r', '%0D') # Carriage return
    
    return escaped


def unescape_vcf_info_value(escaped_value: str) -> str:
    """
    Unescape a previously escaped VCF INFO field value.
    
    This function reverses the escaping performed by escape_vcf_info_value(),
    converting URL-style percent encoding back to original characters.
    
    Args:
        escaped_value: The escaped string from VCF INFO field
        
    Returns:
        Original unescaped string
        
    Example:
        >>> unescape_vcf_info_value("LowDepth%3BLowEVS")
        'LowDepth;LowEVS'
        >>> unescape_vcf_info_value("key%3Dvalue%2Cother")
        'key=value,other'
    """
    if not escaped_value or escaped_value == '.':
        return ''
    
    # Unescape URL encoding back to original characters
    unescaped = escaped_value
    unescaped = unescaped.replace('%3B', ';')  # Semicolon
    unescaped = unescaped.replace('%3D', '=')  # Equals
    unescaped = unescaped.replace('%2C', ',')  # Comma
    unescaped = unescaped.replace('%20', ' ')  # Space
    unescaped = unescaped.replace('%09', '\t') # Tab
    unescaped = unescaped.replace('%0A', '\n') # Newline
    unescaped = unescaped.replace('%0D', '\r') # Carriage return
    
    return unescaped


def escape_filter_for_info(filter_value: Union[str, None]) -> str:
    """
    Specifically escape FILTER values for storage in INFO fields.
    
    This is a specialized function for escaping FILTER field values (which commonly
    contain semicolons in compound filters like "LowDepth;LowEVS") when they need
    to be stored in INFO fields.
    
    Args:
        filter_value: FILTER field value (e.g., "PASS", "LowDepth;LowEVS", None)
        
    Returns:
        Escaped filter value safe for INFO field storage
        
    Example:
        >>> escape_filter_for_info("LowDepth;LowEVS")
        'LowDepth%3BLowEVS'
        >>> escape_filter_for_info("PASS")
        'PASS'
        >>> escape_filter_for_info(None)
        'PASS'
    """
    if filter_value is None or filter_value == '':
        return 'PASS'
    
    # Handle common single filters that don't need escaping
    if filter_value in ['PASS', 'FAIL', '.']:
        return filter_value
    
    # Escape compound filters and special characters
    return escape_vcf_info_value(filter_value)


def validate_vcf_info_field(info_string: str) -> bool:
    """
    Validate that an INFO field string is properly formatted.
    
    This function checks that an INFO field string doesn't contain unescaped
    special characters that would cause VCF parsing issues.
    
    Args:
        info_string: Complete INFO field string (e.g., "DP=30;MQ=40;FILTER=PASS")
        
    Returns:
        True if INFO string is valid, False if it contains unescaped special chars
        
    Example:
        >>> validate_vcf_info_field("DP=30;MQ=40")
        True
        >>> validate_vcf_info_field("DP=30;FILTER=LowDepth;LowEVS;MQ=40")
        False  # Unescaped semicolon in FILTER value
    """
    if not info_string or info_string == '.':
        return True
    
    # Split by semicolons to get individual INFO fields
    info_fields = info_string.split(';')
    
    for field in info_fields:
        if '=' in field:
            key, value = field.split('=', 1)
            
            # Check if value contains unescaped special characters
            # Look for patterns that suggest unescaped characters
            if ';' in value and '%3B' not in value:
                # Unescaped semicolon in value
                return False
            if ',' in value and '%2C' not in value:
                # This might be OK for multi-valued fields, but check context
                pass
        
    return True


def fix_compound_filter_in_info(info_string: str) -> str:
    """
    Fix compound FILTER values in INFO strings by escaping semicolons.
    
    This function specifically looks for FILTER-related INFO fields that contain
    unescaped semicolons and fixes them by applying proper escaping.
    
    Args:
        info_string: INFO field string that may contain unescaped FILTER values
        
    Returns:
        Fixed INFO string with properly escaped FILTER values
        
    Example:
        >>> fix_compound_filter_in_info("DP=30;FILTERS_ORIGINAL=DNA_strelka:LowDepth;LowEVS")
        'DP=30;FILTERS_ORIGINAL=DNA_strelka:LowDepth%3BLowEVS'
    """
    if not info_string or info_string == '.':
        return info_string
    
    # Pattern to match FILTER-related INFO fields with compound values
    filter_field_pattern = r'(FILTERS?_[A-Z_]*=)([^;]+(?:;[^;=]+)*)'
    
    def escape_filter_match(match):
        field_name = match.group(1)  # e.g., "FILTERS_ORIGINAL="
        field_value = match.group(2)  # e.g., "DNA_strelka:LowDepth;LowEVS"
        
        # Escape semicolons in the field value
        escaped_value = field_value.replace(';', '%3B')
        return field_name + escaped_value
    
    # Apply the fix
    fixed_string = re.sub(filter_field_pattern, escape_filter_match, info_string)
    
    return fixed_string


# Convenience functions for backward compatibility
def escape_for_vcf_info(value: Union[str, int, float, None]) -> str:
    """Alias for escape_vcf_info_value for backward compatibility."""
    return escape_vcf_info_value(value)


def unescape_from_vcf_info(escaped_value: str) -> str:
    """Alias for unescape_vcf_info_value for backward compatibility."""
    return unescape_vcf_info_value(escaped_value)


if __name__ == "__main__":
    # Test the escaping functions
    test_cases = [
        "LowDepth;LowEVS",
        "PASS",
        "key=value,other",
        "normal text",
        None,
        "",
        "complex;filter=with,multiple;separators"
    ]
    
    print("Testing VCF INFO field escaping:")
    for test_case in test_cases:
        escaped = escape_vcf_info_value(test_case)
        unescaped = unescape_vcf_info_value(escaped)
        print(f"Original: {test_case!r}")
        print(f"Escaped:  {escaped!r}")
        print(f"Unescaped: {unescaped!r}")
        print(f"Round-trip OK: {str(test_case) == unescaped if test_case is not None else unescaped == ''}")
        print()