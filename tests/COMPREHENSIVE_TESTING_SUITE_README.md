# Comprehensive Testing Suite for VCF Realignment Optimization

## Overview

This comprehensive testing suite implements property-based testing and debugging tools for the VCF Realignment Optimization feature. The suite addresses all requirements from the design document and provides robust validation of the system's correctness properties.

## Components

### 1. Property-Based Test Framework

**Location:** `tests/property_test_framework.nf`

**Purpose:** Implements property-based testing with minimum 100 iterations per property to validate universal correctness properties.

**Key Features:**
- Test data generators for CRAM, VCF, and metadata
- Configurable iteration counts (minimum 100)
- Reproducible test data with seeded random generation
- Comprehensive property validation

**Properties Tested:**
1. **StackOverflow Prevention** - Validates Requirements 1.1, 1.2, 1.3
2. **Null Condition Handling** - Validates Requirements 1.5
3. **Channel Join Integrity** - Validates Requirements 3.1, 3.2
4. **File Path Validation** - Validates Requirements 3.3
5. **Command Construction Validation** - Validates Requirements 5.5

### 2. Individual Property Test Files

**Location:** `tests/property_tests/`

Each property has its own dedicated test file with multiple test scenarios:

- `stackoverflow_prevention.nf.test` - Tests for StackOverflow prevention
- `null_condition_handling.nf.test` - Tests for null condition handling
- `channel_join_integrity.nf.test` - Tests for channel join operations
- `file_path_validation.nf.test` - Tests for file validation
- `command_construction.nf.test` - Tests for command validation

**Test Configuration:**
- Minimum 100 iterations per property test
- Tagged with format: **Feature: vcf-realignment-optimization, Property {number}: {property_text}**
- Configurable test parameters (iterations, timeouts, resources)

### 3. Debug and Diagnostic Tools

**Location:** `tests/debug_tools/`

**Components:**
- `channel_inspector.nf` - Channel content and flow analysis
- `meta_analyzer.nf` - Metadata structure analysis and circular reference detection
- `performance_profiler.nf` - Resource usage and performance monitoring
- `debug_toolkit.nf` - Comprehensive debugging workflow

**Capabilities:**
- **Channel Inspection:** Analyze channel contents, detect circular references, monitor memory usage
- **Meta Analysis:** Deep metadata structure analysis, validation, sanitization testing
- **Performance Profiling:** Workflow execution profiling, throughput monitoring, bottleneck detection
- **Integrated Debugging:** Comprehensive debugging combining all tools

### 4. Test Runners and Configuration

**Property Test Runner:** `tests/run_property_tests.sh`
- Automated execution of all property tests
- Configurable iterations and timeouts
- Comprehensive reporting
- Test environment setup and validation

**Debug Tools Runner:** `tests/run_debug_tools.sh`
- Individual tool execution
- Comprehensive debugging workflows
- Test data generation
- Report generation

**Configuration:** `tests/config/property_test.config`
- Property test parameters
- Resource limits
- Test data paths
- Logging configuration

## Usage

### Running Property Tests

```bash
# Setup test environment
./tests/run_property_tests.sh setup

# Run all property tests
./tests/run_property_tests.sh run

# Generate test report
./tests/run_property_tests.sh report
```

### Running Debug Tools

```bash
# Run all debug tools
./tests/run_debug_tools.sh all

# Run specific tools
./tests/run_debug_tools.sh channel_inspector test_name /path/to/vcf /path/to/cram
./tests/run_debug_tools.sh meta_analyzer test_name /path/to/input
./tests/run_debug_tools.sh performance_profiler test_name /path/to/input 100
./tests/run_debug_tools.sh debug_toolkit test_name /path/to/vcf /path/to/cram

# Generate debug report
./tests/run_debug_tools.sh report
```

### Running Individual Property Tests

```bash
# Run specific property test
nf-test test tests/property_tests/stackoverflow_prevention.nf.test

# Run with custom configuration
nf-test test tests/property_tests/channel_join_integrity.nf.test --config tests/config/property_test.config
```

## Test Data Generation

The framework includes sophisticated test data generators:

### Metadata Generator
- Generates random metadata with controlled patient ID distribution
- Supports various sample statuses (DNA normal, DNA tumor, RNA)
- Creates realistic metadata structures for testing

### VCF Generator
- Generates valid VCF content with random variants
- Configurable sample names and variant counts
- Proper VCF header formatting

### CRAM Generator
- Simulates CRAM file metadata
- Configurable file sizes and read counts
- Mock file creation for testing

## Property Test Configuration

### Default Parameters
- **Iterations:** 100 minimum per property
- **Timeout:** 30 minutes per test
- **Max Samples:** 10 per test
- **Max Patients:** 5 per test
- **Resources:** 2 CPUs, 4GB memory

### Customization
All parameters can be customized through the configuration file or command-line parameters.

## Requirements Validation

### Property-Based Testing Requirements
- ✅ **All property requirements** - Comprehensive property test framework
- ✅ **Minimum 100 iterations** - Configurable with 100 as minimum
- ✅ **Test data generators** - CRAM, VCF, and metadata generators
- ✅ **nf-test configuration** - Proper nf-test setup and configuration

### Debugging and Diagnostic Requirements
- ✅ **Requirements 6.5** - Channel inspection utilities
- ✅ **Requirements 6.5** - Meta data structure analysis tools
- ✅ **Requirements 6.5** - Performance profiling capabilities

### Correctness Properties Validation
- ✅ **Property 1:** StackOverflow Prevention (Requirements 1.1, 1.2, 1.3)
- ✅ **Property 2:** Null Condition Handling (Requirements 1.5)
- ✅ **Property 3:** Channel Join Integrity (Requirements 3.1, 3.2)
- ✅ **Property 4:** File Path Validation (Requirements 3.3)
- ✅ **Property 7:** Command Construction Validation (Requirements 5.5)

## Integration with Existing Tests

The comprehensive testing suite integrates with existing test infrastructure:

- Uses existing `nf-test.config` configuration
- Leverages existing test data paths
- Compatible with existing CI/CD pipelines
- Extends existing integration tests

## Reporting and Analysis

### Property Test Reports
- Test execution summaries
- Property validation results
- Performance metrics
- Failure analysis

### Debug Tool Reports
- Channel analysis results
- Metadata health assessments
- Performance profiles
- Bottleneck identification

### Comprehensive Reports
- Combined property and debug analysis
- Requirements traceability
- Issue identification and recommendations
- Overall system health assessment

## Maintenance and Extension

### Adding New Properties
1. Define property in design document
2. Add property test to framework
3. Create dedicated test file
4. Update test runner
5. Update documentation

### Adding New Debug Tools
1. Create tool in `debug_tools/` directory
2. Integrate with debug toolkit
3. Update debug runner
4. Add usage documentation

### Customizing Test Data
1. Modify generators in `property_test_framework.nf`
2. Update configuration parameters
3. Add new data types as needed

## Best Practices

### Property Test Design
- Focus on universal properties, not specific examples
- Use realistic test data generation
- Validate essential system behaviors
- Include edge case handling

### Debug Tool Usage
- Use channel inspector for flow analysis
- Use meta analyzer for circular reference detection
- Use performance profiler for bottleneck identification
- Use comprehensive toolkit for full analysis

### Test Maintenance
- Regular execution of property tests
- Monitor test performance and adjust iterations
- Update test data generators as system evolves
- Maintain requirements traceability

## Conclusion

This comprehensive testing suite provides robust validation of the VCF Realignment Optimization feature through property-based testing and comprehensive debugging tools. It addresses all requirements from the design document and provides a solid foundation for ensuring system correctness and maintainability.