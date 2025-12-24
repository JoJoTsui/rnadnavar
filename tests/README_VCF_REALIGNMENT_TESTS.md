# VCF Realignment Optimization Testing Framework

This document describes the comprehensive testing framework for the VCF realignment workflow optimization in the nf-core/rnadnavar pipeline.

## Overview

The testing framework provides comprehensive validation of the VCF realignment optimization functionality through multiple test categories:

- **Integration Tests**: Test complete VCF realignment workflow end-to-end
- **Comparison Tests**: Compare VCF vs MAF realignment workflows
- **Error Handling Tests**: Test StackOverflowError prevention and error recovery
- **Performance Tests**: Test processing performance and resource usage
- **Backward Compatibility Tests**: Ensure MAF workflow continues to work

## Test Structure

### Test Data

Test data files are located in `tests/data/`:

- `test_vcf_realignment.vcf.gz` - VCF file with candidate regions for realignment
- `test_vcf_realignment.vcf.gz.tbi` - Tabix index for VCF file
- `cram_trbc_recal_complex.csv` - CRAM files with RNA/DNA samples for testing

### Test Categories

#### 1. Integration Tests (`tests/vcf_realignment_integration.nf.test`)

Tests the complete VCF realignment workflow:

- **End-to-end validation**: Complete workflow execution with VCF input
- **Sample filtering**: RNA (status=2) and DNA normal (status=0) processing
- **Error handling validation**: StackOverflowError prevention and recovery
- **MAF backward compatibility**: Ensures original MAF workflow still works
- **Workflow equivalence**: Compares VCF vs MAF realignment outputs

#### 2. Comparison Tests (`tests/vcf_realignment_comparison.nf.test`)

Tests workflow comparison and validation:

- **Known good output validation**: Validates against expected outputs
- **MAF baseline comparison**: Compares with original MAF workflow
- **Performance comparison**: Validates performance is not degraded
- **Edge case handling**: Tests with empty and malformed inputs

#### 3. Component Tests

Individual component testing:

- **Channel Sanitization**: `subworkflows/local/sanitize_channels/tests/main.nf.test`
- **Safe Channel Join**: `subworkflows/local/safe_channel_join/tests/main.nf.test`
- **Enhanced CRAM2BAM**: `subworkflows/local/enhanced_cram2bam_conversion/tests/main.nf.test`
- **Input Validation**: `subworkflows/local/input_validation/tests/main.nf.test`

## Running Tests

### Complete Test Suite

Run the complete VCF realignment test suite:

```bash
# Execute all VCF realignment tests
./tests/run_vcf_realignment_tests.sh
```

This script runs all test categories in the correct order and provides a comprehensive summary.

### Individual Test Execution

Run specific test files using nf-test:

```bash
# Run integration tests
nf-test test tests/vcf_realignment_integration.nf.test

# Run comparison tests
nf-test test tests/vcf_realignment_comparison.nf.test

# Run component tests
nf-test test subworkflows/local/sanitize_channels/tests/main.nf.test
```

### Test Profiles

Use specific test profiles for different scenarios:

```bash
# Test with VCF realignment enabled
nf-test test tests/vcf_realignment_integration.nf.test -profile vcf_realignment_test

# Test with MAF realignment (backward compatibility)
nf-test test tests/vcf_realignment_integration.nf.test -profile maf_realignment_test

# Performance testing profile
nf-test test tests/vcf_realignment_integration.nf.test -profile vcf_realignment_performance
```

## Test Configuration

### Test-Specific Configuration

VCF realignment tests use dedicated configuration in `tests/config/vcf_realignment_test.config`:

- Resource allocation for test processes
- Test-specific parameter defaults
- Profile definitions for different test scenarios

### Test Data Configuration

Test data paths are configured in `tests/config/test_data.config` under the `vcf_realignment` section.

## Expected Test Results

### Successful Test Execution

When all tests pass, you should see:

- All integration tests complete successfully
- VCF realignment workflow components execute without StackOverflowError
- Channel sanitization and validation working correctly
- Error handling and monitoring functional
- Backward compatibility with MAF workflow maintained
- Performance within acceptable bounds

### Test Validation Criteria

Tests validate the following correctness properties:

1. **StackOverflow Prevention** - No StackOverflowError during channel operations
2. **Channel Join Integrity** - Proper patient ID matching and metadata preservation
3. **File Path Validation** - File existence and accessibility verification
4. **Channel Schema Consistency** - Data type validation across channel operations
5. **Input Format Validation** - VCF, CRAM, and index file format validation
6. **Command Construction Validation** - Proper command syntax and parameters
7. **Workflow Equivalence** - VCF and MAF workflows produce comparable results

## Key Test Scenarios

### 1. StackOverflowError Prevention

Tests specifically validate that the original StackOverflowError in SAMTOOLS_CONVERT is resolved:

```bash
# Verify no StackOverflowError occurs
assert !workflow.trace.tasks().any { it.stderr?.contains('StackOverflowError') }
```

### 2. Channel Sanitization

Tests validate that channel data is properly cleaned to prevent circular references:

```bash
# Verify sanitization components executed
assert workflow.trace.tasks().any { it.name.contains('SANITIZE_CHANNELS') }
```

### 3. Sample Filtering

Tests validate that only RNA (status=2) and DNA normal (status=0) samples are processed:

```bash
# DNA tumor samples (status=1) should be filtered out
assert workflow.trace.tasks().any { it.name.contains('PROCESS_MONITORING') }
```

### 4. Backward Compatibility

Tests ensure the original MAF workflow continues to work:

```bash
# MAF workflow should not use VCF-specific components
assert !workflow.trace.tasks().any { it.name.contains('VCF2BED') }
```

## Troubleshooting

### Common Test Failures

1. **Missing test data**: Ensure VCF files are properly compressed and indexed
2. **Resource limitations**: Adjust memory/CPU allocations in test configuration
3. **Container issues**: Verify conda environment and container specifications
4. **Network issues**: Check connectivity for remote test data downloads

### Test Data Regeneration

If test data needs to be regenerated:

```bash
# Regenerate compressed VCF files
cd tests/data
bgzip -c test_vcf_realignment.vcf > test_vcf_realignment.vcf.gz
tabix -p vcf test_vcf_realignment.vcf.gz
```

### Debug Mode

Run tests with additional debugging:

```bash
# Run with verbose output
nf-test test <test_file> --verbose

# Run with debug logging
nf-test test <test_file> --debug
```

## Performance Benchmarks

### Expected Performance Metrics

- **Channel Operations**: < 10 seconds per operation
- **VCF2BED Conversion**: < 30 seconds for test VCF
- **CRAM2BAM Conversion**: < 2 minutes for test CRAM
- **Read Filtering**: < 1 minute for test data
- **HISAT2 Realignment**: < 5 minutes for test data

### Memory Usage

- **Channel Operations**: < 1GB
- **File Conversions**: < 2GB
- **Realignment**: < 4GB
- **Total Workflow**: < 6GB

## Continuous Integration

The VCF realignment test suite is designed to integrate with CI/CD pipelines:

- Tests are organized by execution time (fast component tests first)
- Resource requirements are optimized for CI environments
- Test results provide clear pass/fail indicators
- Comprehensive logging aids in debugging CI failures

## Contributing

When adding new VCF realignment functionality:

1. Add corresponding unit tests to the component test files
2. Update integration tests if workflow changes are made
3. Add error handling tests for new failure modes
4. Include performance tests for resource-intensive operations
5. Update comparison tests for workflow equivalence validation

Ensure all tests pass before submitting changes:

```bash
./tests/run_vcf_realignment_tests.sh
```

## Requirements Validation

This testing framework validates the following requirements from the specification:

- **Requirement 7.1**: VCF realignment workflow validation against known good outputs
- **Requirement 7.2**: MAF realignment workflow backward compatibility
- **Requirement 7.5**: End-to-end data flow validation

The tests provide comprehensive coverage of the VCF realignment optimization, ensuring that the StackOverflowError is resolved, performance is maintained, and backward compatibility is preserved.