# Test inputs

These files come from the Benchling-Open-Source/allotropy test data
(https://github.com/Benchling-Open-Source/allotropy), MIT License,
Copyright (c) Benchling. Original paths are relative to the allotropy
repository root.

- `openlab-sirius-single-injection.rslt`: `tests/parsers/agilent_openlab_cds/testdata/Sirius-single-injection-2023-09-01 07-52-44-04-00.rslt`, unchanged.
- `openlab-sirius.rslt`: `tests/parsers/agilent_openlab_cds/testdata/Sirius-2023-09-01 07-52-44-04-00.rslt`, unchanged.
- `unicorn-single-uv.zip`: `tests/parsers/cytiva_unicorn/testdata/unicorn_single_uv.zip`, unchanged.
- `unicorn-multiple-chromatograms.zip`: `tests/parsers/cytiva_unicorn/testdata/unicorn_multiple_chromatograms.zip`, unchanged.
- `unicorn-1-trimmed.zip`: `tests/parsers/cytiva_unicorn/testdata/unicorn_1.zip`, only the members ending in Result.xml, Chrom.1.Xml, SystemData.zip, InstrumentConfigurationData.zip, EvaluationLog.xml (and the folder entry) kept, contents and compression unchanged.
- `empower-example-01.head.json`: `tests/parsers/benchling_empower/testdata/input/example_01.json`, first 32768 bytes.
- `empower-blanks-and-stds.json`: `tests/parsers/benchling_empower/testdata/input/example_with_blanks_and_stds.json`, re-serialized without indentation (same JSON value).
- `chromeleon-multi-signal.json`: `tests/parsers/benchling_chromeleon/testdata/input/chromeleon_example_multi_signal.json`, unchanged.
- `chromeleon-null-device-info.json`: `tests/parsers/benchling_chromeleon/testdata/input/chromeleon_example_null_device_info.json`, unchanged.
- `chromeleon-example.head.json`: `tests/parsers/benchling_chromeleon/testdata/input/chromeleon_example.json`, first 32768 bytes.
- `tapestation-example-01.xml`: `tests/parsers/agilent_tapestation_analysis/testdata/agilent_tapestation_analysis_example_01.xml`, unchanged.
- `tapestation-example-03.xml`: `tests/parsers/agilent_tapestation_analysis/testdata/agilent_tapestation_analysis_example_03.xml`, unchanged.
- `tapestation-example-02.head.xml`: `tests/parsers/agilent_tapestation_analysis/testdata/agilent_tapestation_analysis_example_02_error.xml`, first 32768 bytes.
- `insight-transposed-flowcells.xlsx`: `tests/parsers/cytiva_biacore_insight/testdata/Biacore_Insight_Transposed_FlowCells.xlsx`, unchanged.
- `t200-control-immobilization.blr.sparse.json`: `tests/parsers/cytiva_biacore_t200_control/testdata/ED_Fig.6a_immobilization Her2-Her3.blr`, sparse recording: only the byte ranges the sniffer reads, as hex, with the file size.
