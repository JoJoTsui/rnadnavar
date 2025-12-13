#!/usr/bin/env python3
"""
REDIportal Field Mapping Module

This module provides accurate field mapping functionality for converting REDIportal
database fields to VCF INFO fields with proper escaping and validation.

Key Features:
- Static methods for mapping each REDIportal field to VCF INFO field
- Proper escaping of special characters according to VCF specification
- Validation of field values
- Consistent handling of empty/missing values

Author: RNA Editing Enhancement Pipeline
Date: 2025-12-13
"""

import logging
import re
from typing import Optional

logger = logging.getLogger(__name__)

class REDIportalFieldMapper:
    """
    Static class for mapping REDIportal fields to VCF INFO fields.
    
    This class provides static methods for converting each REDIportal field
    to its corresponding VCF INFO field with proper escaping and validation.
    """
    
    # VCF INFO field names for REDIportal data
    INFO_FIELD_NAMES = {
        'ACCESSION': 'REDI_ACCESSION',
        'DB': 'REDI_DB',
        'TYPE': 'REDI_TYPE',
        'REPEAT': 'REDI_REPEAT',
        'STRAND': 'REDI_STRAND'
    }
    
    # VCF header definitions for REDIportal INFO fields
    INFO_FIELD_DEFINITIONS = {
        'REDI_ACCESSION': '##INFO=<ID=REDI_ACCESSION,Number=1,Type=String,Description="REDIportal accession identifier">',
        'REDI_DB': '##INFO=<ID=REDI_DB,Number=1,Type=String,Description="REDIportal database source">',
        'REDI_TYPE': '##INFO=<ID=REDI_TYPE,Number=1,Type=String,Description="REDIportal editing type classification">',
        'REDI_REPEAT': '##INFO=<ID=REDI_REPEAT,Number=1,Type=String,Description="REDIportal repeat element annotation">',
        'REDI_STRAND': '##INFO=<ID=REDI_STRAND,Number=1,Type=String,Description="REDIportal strand information">'
    }
    
    @staticmethod
    def escape_info_value(value: str) -> str:
        """
        Escape special characters in INFO field values according to VCF specification.
        
        VCF specification requires escaping of certain characters in INFO field values:
        - Semicolon (;) is the field separator
        - Equals (=) is the key-value separator  
        - Comma (,) is the value separator for multi-value fields
        - Whitespace characters should be replaced with underscores
        
        Args:
            value: Raw value to escape
            
        Returns:
            Escaped value safe for VCF INFO field
        """
        if not value or value.strip() == '' or value == '.':
            return '.'
        
        # Strip whitespace first
        value = value.strip()
        
        # Replace VCF special characters with URL encoding
        value = value.replace(';', '%3B')  # Semicolon is field separator
        value = value.replace('=', '%3D')  # Equals is key-value separator
        value = value.replace(',', '%2C')  # Comma is value separator
        
        # Replace whitespace characters with underscores
        value = value.replace(' ', '_')    # Spaces to underscores
        value = value.replace('\t', '_')   # Tabs to underscores
        value = value.replace('\n', '_')   # Newlines to underscores
        value = value.replace('\r', '_')   # Carriage returns to underscores
        
        # Replace other problematic characters
        value = value.replace('"', '%22')  # Double quotes
        value = value.replace("'", '%27')  # Single quotes
        value = value.replace('\\', '%5C') # Backslashes
        
        # Remove any remaining control characters
        value = re.sub(r'[\x00-\x1f\x7f-\x9f]', '_', value)
        
        return value
    
    @staticmethod
    def validate_field_value(value: str, field_name: str) -> bool:
        """
        Validate that a field value is appropriate for VCF INFO field.
        
        Args:
            value: Field value to validate
            field_name: Name of the field for logging
            
        Returns:
            True if value is valid, False otherwise
        """
        if not value or value.strip() == '':
            return True  # Empty values are valid (will be converted to '.')
        
        # Check for extremely long values that might cause issues
        if len(value) > 1000:
            logger.warning(f"Field {field_name} has very long value ({len(value)} chars), truncating")
            return False
        
        # Check for binary data or non-printable characters
        try:
            value.encode('utf-8')
        except UnicodeEncodeError:
            logger.warning(f"Field {field_name} contains non-UTF8 characters")
            return False
        
        return True
    
    @staticmethod
    def map_accession_field(accession: str) -> str:
        """
        Map REDIportal Accession field to REDI_ACCESSION INFO field.
        
        The Accession field contains REDIportal database identifiers that need
        proper escaping for VCF INFO field usage.
        
        Args:
            accession: Raw Accession field value from REDIportal
            
        Returns:
            Escaped value suitable for REDI_ACCESSION INFO field
        """
        if not REDIportalFieldMapper.validate_field_value(accession, 'Accession'):
            return '.'
        
        escaped_value = REDIportalFieldMapper.escape_info_value(accession)
        
        # Log debug info for non-trivial mappings
        if escaped_value != accession and escaped_value != '.':
            logger.debug(f"Accession field escaped: '{accession}' -> '{escaped_value}'")
        
        return escaped_value
    
    @staticmethod
    def map_db_field(db: str) -> str:
        """
        Map REDIportal db field to REDI_DB INFO field.
        
        The db field contains database source information that needs
        proper escaping for VCF INFO field usage.
        
        Args:
            db: Raw db field value from REDIportal
            
        Returns:
            Escaped value suitable for REDI_DB INFO field
        """
        if not REDIportalFieldMapper.validate_field_value(db, 'db'):
            return '.'
        
        escaped_value = REDIportalFieldMapper.escape_info_value(db)
        
        # Log debug info for non-trivial mappings
        if escaped_value != db and escaped_value != '.':
            logger.debug(f"db field escaped: '{db}' -> '{escaped_value}'")
        
        return escaped_value
    
    @staticmethod
    def map_type_field(type_val: str) -> str:
        """
        Map REDIportal type field to REDI_TYPE INFO field.
        
        The type field contains RNA editing type classification that needs
        proper escaping for VCF INFO field usage.
        
        Args:
            type_val: Raw type field value from REDIportal
            
        Returns:
            Escaped value suitable for REDI_TYPE INFO field
        """
        if not REDIportalFieldMapper.validate_field_value(type_val, 'type'):
            return '.'
        
        escaped_value = REDIportalFieldMapper.escape_info_value(type_val)
        
        # Log debug info for non-trivial mappings
        if escaped_value != type_val and escaped_value != '.':
            logger.debug(f"type field escaped: '{type_val}' -> '{escaped_value}'")
        
        return escaped_value
    
    @staticmethod
    def map_repeat_field(repeat: str) -> str:
        """
        Map REDIportal repeat field to REDI_REPEAT INFO field.
        
        The repeat field contains repeat element annotation that needs
        proper escaping for VCF INFO field usage.
        
        Args:
            repeat: Raw repeat field value from REDIportal
            
        Returns:
            Escaped value suitable for REDI_REPEAT INFO field
        """
        if not REDIportalFieldMapper.validate_field_value(repeat, 'repeat'):
            return '.'
        
        escaped_value = REDIportalFieldMapper.escape_info_value(repeat)
        
        # Log debug info for non-trivial mappings
        if escaped_value != repeat and escaped_value != '.':
            logger.debug(f"repeat field escaped: '{repeat}' -> '{escaped_value}'")
        
        return escaped_value
    
    @staticmethod
    def map_strand_field(strand: str) -> str:
        """
        Map REDIportal Strand field to REDI_STRAND INFO field.
        
        The Strand field contains strand information that needs
        proper escaping for VCF INFO field usage.
        
        Args:
            strand: Raw Strand field value from REDIportal
            
        Returns:
            Escaped value suitable for REDI_STRAND INFO field
        """
        if not REDIportalFieldMapper.validate_field_value(strand, 'Strand'):
            return '.'
        
        escaped_value = REDIportalFieldMapper.escape_info_value(strand)
        
        # Log debug info for non-trivial mappings
        if escaped_value != strand and escaped_value != '.':
            logger.debug(f"Strand field escaped: '{strand}' -> '{escaped_value}'")
        
        return escaped_value
    
    @staticmethod
    def get_info_field_definitions() -> dict:
        """
        Get VCF header definitions for all REDIportal INFO fields.
        
        Returns:
            Dictionary mapping INFO field names to their VCF header definitions
        """
        return REDIportalFieldMapper.INFO_FIELD_DEFINITIONS.copy()
    
    @staticmethod
    def get_info_field_names() -> dict:
        """
        Get mapping of REDIportal field types to VCF INFO field names.
        
        Returns:
            Dictionary mapping field types to INFO field names
        """
        return REDIportalFieldMapper.INFO_FIELD_NAMES.copy()
    
    @staticmethod
    def map_all_fields(accession: str, db: str, type_val: str, repeat: str, strand: str) -> dict:
        """
        Map all REDIportal fields to their corresponding VCF INFO fields.
        
        This is a convenience method that maps all fields at once and returns
        a dictionary suitable for building VCF INFO field strings.
        
        Args:
            accession: Raw Accession field value
            db: Raw db field value
            type_val: Raw type field value
            repeat: Raw repeat field value
            strand: Raw Strand field value
            
        Returns:
            Dictionary mapping INFO field names to escaped values
        """
        return {
            'REDI_ACCESSION': REDIportalFieldMapper.map_accession_field(accession),
            'REDI_DB': REDIportalFieldMapper.map_db_field(db),
            'REDI_TYPE': REDIportalFieldMapper.map_type_field(type_val),
            'REDI_REPEAT': REDIportalFieldMapper.map_repeat_field(repeat),
            'REDI_STRAND': REDIportalFieldMapper.map_strand_field(strand)
        }
    
    @staticmethod
    def build_info_string(field_mappings: dict) -> str:
        """
        Build VCF INFO field string from field mappings.
        
        Args:
            field_mappings: Dictionary mapping INFO field names to values
            
        Returns:
            VCF INFO field string (semicolon-separated key=value pairs)
        """
        info_parts = []
        
        for field_name, value in field_mappings.items():
            if value and value != '.':
                info_parts.append(f"{field_name}={value}")
        
        return ';'.join(info_parts) if info_parts else '.'


# Convenience functions for backward compatibility
def map_accession_field(accession: str) -> str:
    """Convenience function for mapping Accession field."""
    return REDIportalFieldMapper.map_accession_field(accession)

def map_db_field(db: str) -> str:
    """Convenience function for mapping db field."""
    return REDIportalFieldMapper.map_db_field(db)

def map_type_field(type_val: str) -> str:
    """Convenience function for mapping type field."""
    return REDIportalFieldMapper.map_type_field(type_val)

def map_repeat_field(repeat: str) -> str:
    """Convenience function for mapping repeat field."""
    return REDIportalFieldMapper.map_repeat_field(repeat)

def map_strand_field(strand: str) -> str:
    """Convenience function for mapping Strand field."""
    return REDIportalFieldMapper.map_strand_field(strand)

def escape_info_value(value: str) -> str:
    """Convenience function for escaping INFO field values."""
    return REDIportalFieldMapper.escape_info_value(value)