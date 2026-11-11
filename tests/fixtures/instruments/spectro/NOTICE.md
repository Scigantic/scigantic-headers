# Spectrophotometer export fixtures

These files come from the tests of Benchling-Open-Source/allotropy
(https://github.com/Benchling-Open-Source/allotropy), MIT License,
Copyright (c) Benchling. They are copied unmodified from commit
88f59dc1f0bf8efbd01cf18ae3140cf4cd1aeb60, renamed as listed below.

| file | original path in allotropy |
| --- | --- |
| visionlite_fixed_result_multi.csv | tests/parsers/thermo_fisher_visionlite/testdata/Fixed_with_Result_Multi.csv |
| visionlite_quant.csv | tests/parsers/thermo_fisher_visionlite/testdata/Thermo_VISIONlite_example_quant.csv |
| visionlite_scan.csv | tests/parsers/thermo_fisher_visionlite/testdata/Thermo_VISIONlite_example_scan.csv |
| visionlite_rate.csv | tests/parsers/thermo_fisher_visionlite/testdata/exclude/Thermo_VISIONlite_example_rate.csv |
| visionlite_fixed_tab_utf16.csv | tests/parsers/thermo_fisher_visionlite/testdata/Fixed_tab_sep.csv |
| genesys_on_board_fixed.csv | tests/parsers/thermo_fisher_genesys_on_board/testdata/Genesys150_ExampleFile_Fixed_10_2_2024_1_57_42_PM.csv |
| nanodrop_8000_protein.txt | tests/parsers/thermo_fisher_nanodrop_8000/testdata/Thermo_NanoDrop_8000_example01.txt |
| nanodrop_8000_dna.txt | tests/parsers/thermo_fisher_nanodrop_8000/testdata/Thermo_NanoDrop_8000_example03.txt |
| nanodrop_8000_mrna.txt | tests/parsers/thermo_fisher_nanodrop_8000/testdata/Thermo_NanoDrop_mRNA.txt |
| nanodrop_eight_dsdna.txt | tests/parsers/thermo_fisher_nanodrop_eight/testdata/thermo_nanodrop_eight_example01.txt |
| nanodrop_eight_rna.txt | tests/parsers/thermo_fisher_nanodrop_eight/testdata/thermo_nanodrop_eight_RNA.txt |
| nanodrop_eight_tab_aligned.txt | tests/parsers/thermo_fisher_nanodrop_eight/testdata/thermo_nanodrop_eight_tab_aligned_header.txt |
| nanodrop_one_dsdna.csv | tests/parsers/thermo_fisher_nanodrop_one/testdata/dsDNA 7_19_2023 3_27_29 PM.csv |
| nanodrop_one_dsdna.xlsx | tests/parsers/thermo_fisher_nanodrop_one/testdata/example1.xlsx |
| qubit4.csv | tests/parsers/thermo_fisher_qubit4/testdata/thermo_fisher_qubit4_example_1.csv |
| qubit4.xlsx | tests/parsers/thermo_fisher_qubit4/testdata/thermo_fisher_qubit4_example_2.xlsx |
| qubit_flex.csv | tests/parsers/thermo_fisher_qubit_flex/testdata/thermo_fisher_qubit_flex_example_01.csv |
| qubit_flex_sep_header.csv | tests/parsers/thermo_fisher_qubit_flex/testdata/thermo_fisher_qubit_flex_example_sep_header.csv |
| qubit_flex.xlsx | tests/parsers/thermo_fisher_qubit_flex/testdata/thermo_fisher_qubit_flex_example_03.xlsx |
| genesys30_scan.tsv | tests/parsers/thermo_fisher_genesys30/testdata/thermo_fisher_genesys30_example_01.tsv |
| genesys30_scan.csv | tests/parsers/thermo_fisher_genesys30/testdata/thermo_fisher_genesys30_example_02.csv |

spectro-cases.json records, for each file, the expected decode and the values
allotropy's expected ASM output (the .json beside the original input) states for
the same file, which the tests cross-check.
