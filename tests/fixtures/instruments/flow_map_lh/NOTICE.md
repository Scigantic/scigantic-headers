# Fixture provenance

These files come from the test data of Benchling-Open-Source/allotropy
(https://github.com/Benchling-Open-Source/allotropy), MIT License,
Copyright (c) 2023 Benchling, Inc. Files marked "first N bytes" are the leading bytes of
the original; the cherry pick report is the original with DETAILS rows 41 to 531
removed so it stays small. Everything else is a byte-for-byte copy.

| fixture | original path in allotropy | copied |
| --- | --- | --- |
| flowjo_example_1.head.wsp | tests/parsers/flowjo/testdata/flowjo_example_1.wsp | first 24576 bytes |
| flowjo_example_3.head.wsp | tests/parsers/flowjo/testdata/flowjo_example_3.wsp | first 4096 bytes |
| facsdiva_example_1.head.xml | tests/parsers/bd_biosciences_facsdiva/testdata/facsdiva_example_1.xml | first 4096 bytes |
| facsdiva_example_2.head.xml | tests/parsers/bd_biosciences_facsdiva/testdata/facsdiva_example_2.xml | first 4096 bytes |
| bio-rad_bio-plex_manager_example_01.head.xml | tests/parsers/biorad_bioplex_manager/testdata/bio-rad_bio-plex_manager_example_01.xml | first 8192 bytes |
| luminex_xPONENT_example02.csv | tests/parsers/luminex_xponent/testdata/luminex_xPONENT_example02.csv | whole file |
| luminex_xPONENT_example02_saved.csv | tests/parsers/luminex_xponent/testdata/luminex_xPONENT_example02_saved.csv | whole file |
| luminex_xPONENT_example03.csv | tests/parsers/luminex_xponent/testdata/luminex_xPONENT_example03.csv | whole file |
| luminex_intelliflex_example_01.csv | tests/parsers/luminex_intelliflex/testdata/luminex_intelliflex_example_01.csv | whole file |
| luminex_intelliflex_single_dataset.csv | tests/parsers/luminex_intelliflex/testdata/luminex_intelliflex_single_dataset.csv | whole file |
| luminex_intelliflex_v2_1.csv | tests/parsers/luminex_intelliflex/testdata/luminex_intelliflex_v2_1.csv | whole file |
| luminex_intelliflex_v2_2.csv | tests/parsers/luminex_intelliflex/testdata/luminex_intelliflex_v2_2.csv | whole file |
| luminex_intelliflex_v2_2_excel.xlsx | tests/parsers/luminex_intelliflex/testdata/luminex_intelliflex_v2_2_excel.xlsx | whole file |
| biomek_example_01.csv | tests/parsers/beckman_coulter_biomek/testdata/example_01.csv | whole file |
| biomek_example_03.log | tests/parsers/beckman_coulter_biomek/testdata/example_03.log | whole file |
| biomek_pipetting_example.log | tests/parsers/beckman_coulter_biomek/testdata/pipetting_example.log | whole file |
| biomek_UnifiedTransfer09-14-2011.log | tests/parsers/beckman_coulter_biomek/testdata/UnifiedTransfer09-14-2011.log | whole file |
| echo_reformat_PartialTransfer.csv | tests/parsers/beckman_echo_plate_reformat/testdata/PartialTransfer.csv | whole file |
| echo_reformat_SurveyVolumeExamples.csv | tests/parsers/beckman_echo_plate_reformat/testdata/SurveyVolumeExamples.csv | whole file |
| echo_reformat_ExceptionsOnly.csv | tests/parsers/beckman_echo_plate_reformat/testdata/ExceptionsOnly.csv | whole file |
| echo_reformat_TransferOnly.csv | tests/parsers/beckman_echo_plate_reformat/testdata/TransferOnly.csv | whole file |
| echo_cherry_pick_CherryPick_Transfer_Complete.trimmed.csv | tests/parsers/beckman_echo_cherry_pick/testdata/CherryPick_Transfer_Complete.csv | lines 1-53 and 545-551 (DETAILS rows 41-531 removed) |

The allotropy license (LICENSE.txt):

MIT License

Copyright (c) 2023 Benchling, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
