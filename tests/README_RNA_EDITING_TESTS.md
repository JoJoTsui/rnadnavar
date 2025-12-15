# RNA Editing Annotation Testing Framework

This document describes the comprehensive testing framework for the RNA editing annotation integration in the nf-core/rnadnavar pipeline.

## Overview

The testing framework provides comprehensive validation of the RNA editing annotation functionality through multiple test categories:

- **Unit Tests**: Test individual modules and components
- **Integration Tests**: Test workflow integration and data flow
- **Error Handling Tests**: Test failure scenarios and error recovery
- **Performance Tests**: Test processing performance with different data sizes
- **Comprehensive Tests**: End-to-end validation of complete functionality

## Test Structure

### Test Data

Test data files are located in `tests/data/`:

- `test_rescue.vcf.gz` - Small rescue VCF for basic testing
- `test_rediportal.vcf.gz` - REDIportal database sample for annotation
- `test_large_rescue.vcf.gz` - Larger VCF for performance testing
- `test_malformed.vcf.gz` - Malformed VCF for error handling tests

All VCF files include corresponding tabix indices (`.tbi` files).

### Test Categories

#### 1. Unit Tests (`modules/local/rna_editing_annotation/tests/main.nf.test`)

Tests the RNA editing annotation module in isolation:

- Basic functionality with standard parameters
- Custom parameter handling (`min_rna_support`)
- Performance with large VCF files
- Error handling with malformed input
- Stub execution testing

#### 2. Integration Tests (`subworkflows/local/vcf_rescue_post_processing/tests/main.nf.test`)

Tests the VCF rescue post-processing subworkflow:

- RNA editing enabled workflow
- RNA editing disabled workflow
- Custom parameter propagation
- Performance with large datasets
- Error handling integration

#### 3. Pipeline Integration Tests (`tests/rna_editing_integration.nf.test`)

Tests complete pipeline integration:

- Full pipeline with RNA editing enabled
- Full pipeline with RNA editing disabled
- Backward compatibility validation
- Parameter propagation through workflow layers

#### 4. Error Handling Tests (`tests/rna_editing_error_handling.nf.test`)

Tests error scenarios and recovery:

- Missing REDIportal database handling
- Invalid parameter handling
- Corrupted VCF input handling
- Resource limitation scenarios

#### 5. Performance Tests (`tests/rna_editing_performance.nf.test`)

Tests processing performance:

- Small VCF baseline performance
- Large VCF processing performance
- Memory usage pattern validation
- Multi-sample scalability testing

#### 6. Comprehensive Tests (`tests/rna_editing_comprehensive.nf.test`)

End-to-end validation tests:

- Full workflow validation with all components
- Parameter validation and configuration
- Modular architecture validation
- Data flow integrity verification

## Running Tests

### Individual Test Execution

Run specific test files using nf-test:

```bash
# Run module unit tests
nf-test test modules/local/rna_editing_annotation/tests/main.nf.test

# Run integration tests
nf-test test subworkflows/local/vcf_rescue_post_processing/tests/main.nf.test

# Run pipeline integration tests
nf-test test tests/rna_editing_integration.nf.test
```

### Complete Test Suite

Run the complete RNA editing test suite:

```bash
# Execute all RNA editing tests
./tests/run_rna_editing_tests.sh
```

This script runs all test categories in the correct order and provides a comprehensive summary.

### Test Profiles

Use specific test profiles for different scenarios:

```bash
# Test with RNA editing enabled
nf-test test tests/rna_editing_integration.nf.test -profile rna_editing_test

# Test with RNA editing disabled
nf-test test tests/rna_editing_integration.nf.test -profile rna_editing_disabled

# Performance testing profile
nf-test test tests/rna_editing_performance.nf.test -profile rna_editing_performance
```

## Test Configuration

### Test-Specific Configuration

RNA editing tests use dedicated configuration in `tests/config/rna_editing_test.config`:

- Resource allocation for test processes
- Test-specific parameter defaults
- Profile definitions for different test scenarios

### Test Data Configuration

Test data paths are configured in `tests/config/test_data.config` under the `rna_editing` section.

## Expected Test Results

### Successful Test Execution

When all tests pass, you should see:

- All module unit tests complete successfully
- Integration tests validate proper data flow
- Error handling tests demonstrate graceful failure recovery
- Performance tests complete within reasonable time/memory bounds
- Comprehensive tests validate end-to-end functionality

### Test Validation Criteria

Tests validate the following correctness properties:

1. **Module wrapper consistency** - RNA editing module correctly wraps the annotation script
2. **Workflow integration preservation** - Existing functionality remains intact
3. **Data flow continuity** - Proper channel handling through rescue → annotation → filtering
4. **Annotation preservation in filtering** - RNA editing annotations preserved during filtering
5. **Modular architecture extensibility** - New processing steps can be inserted
6. **Module reusability** - Module can be used in different workflow contexts
7. **Configuration control** - RNA editing can be enabled/disabled via parameters
8. **Parameter customization** - Custom databases and thresholds are supported
9. **Error isolation** - Failures don't crash the entire workflow
10. **Resource management** - Proper cleanup and resource handling
11. **Comprehensive logging integration** - Processing statistics and metrics are logged

## Troubleshooting

### Common Test Failures

1. **Missing test data**: Ensure test VCF files are properly compressed and indexed
2. **Resource limitations**: Adjust memory/CPU allocations in test configuration
3. **Container issues**: Verify conda environment and container specifications
4. **Path issues**: Check that test data paths are correctly specified

### Test Data Regeneration

If test data needs to be regenerated:

```bash
# Regenerate compressed VCF files
cd tests/data
bgzip -c test_rescue.vcf > test_rescue.vcf.gz
tabix -p vcf test_rescue.vcf.gz

bgzip -c test_rediportal.vcf > test_rediportal.vcf.gz
tabix -p vcf test_rediportal.vcf.gz

bgzip -c test_large_rescue.vcf > test_large_rescue.vcf.gz
tabix -p vcf test_large_rescue.vcf.gz
```

### Debug Mode

Run tests with additional debugging:

```bash
# Run with verbose output
nf-test test <test_file> --verbose

# Run with debug logging
nf-test test <test_file> --debug
```

## Continuous Integration

The RNA editing test suite is designed to integrate with CI/CD pipelines:

- Tests are organized by execution time (fast unit tests first)
- Resource requirements are optimized for CI environments
- Test results provide clear pass/fail indicators
- Comprehensive logging aids in debugging CI failures

## Contributing

When adding new RNA editing functionality:

1. Add corresponding unit tests to the module test file
2. Update integration tests if workflow changes are made
3. Add error handling tests for new failure modes
4. Include performance tests for resource-intensive operations
5. Update comprehensive tests for end-to-end validation

Ensure all tests pass before submitting changes:

```bash
./tests/run_rna_editing_tests.sh
```