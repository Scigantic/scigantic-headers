# Fixture provenance

These files come from the Benchling-Open-Source/allotropy test suite
(https://github.com/Benchling-Open-Source/allotropy), MIT License,
Copyright (c) 2023 Benchling, Inc. They were renamed; the original paths are:

- `gen5-absorbance-pathlength.txt`: `tests/parsers/agilent_gen5/testdata/absorbance/endpoint_pathlength_correct_singleplate.txt`
- `gen5-multiple-read-modes.txt`: `tests/parsers/agilent_gen5/testdata/multi_read_modes/multiple_read_modes.txt`
- `gen5-fluorescence-multiplate.txt`: `tests/parsers/agilent_gen5/testdata/fluorescence/exclude/endpoint_multiplate.txt`
- `gen5-kinetic-fibrillation-head.txt`: `tests/parsers/agilent_gen5/testdata/Synergy instrument datafile (Fibrillation data) - TXT format.txt` (first 23766 bytes, cut at a line end)
- `gen5-image-cell-count.txt`: `tests/parsers/agilent_gen5_image/testdata/Cell_Count_DAPI_GFP.txt`
- `gen5-image-comet-two-reads.txt`: `tests/parsers/agilent_gen5_image/testdata/96-Well Trevigen CometAssayCometChip Imaging and Analysis Sample File 23Nov15.txt`
- `gen5-image-6well-brightfield.txt`: `tests/parsers/agilent_gen5_image/testdata/no_results_agilent_gen5_img_example.txt`
- `softmax-fluorescence-utf16.txt`: `tests/parsers/moldev_softmax_pro/testdata/MD_SMP_fluorescence_endpoint_partial_plate_example01.txt`
- `softmax-absorbance-utf8.txt`: `tests/parsers/moldev_softmax_pro/testdata/abs_endpoint_plates.txt`
- `softmax-luminescence-two-plates.txt`: `tests/parsers/moldev_softmax_pro/testdata/MD_SMP_luminescence_endpoint_example03.txt`
- `softmax-kinetic-spectramax340.txt`: `tests/parsers/moldev_softmax_pro/testdata/spectramax340_kinetic_partial_plate.txt`
- `mars-fluorescence-384.csv`: `tests/parsers/bmg_mars/testdata/16-02-29 14-34-46 Transcreener ADP2 FI.csv`
- `mars-absorbance-384.csv`: `tests/parsers/bmg_mars/testdata/16-03-03 16-54-03 472 ABS 384 QC.csv`
- `mars-luminescence-1536-partial.csv`: `tests/parsers/bmg_mars/testdata/24-10-16 15-09-52 Luminescence 1536.csv`
- `mars-absorbance-96-two-numbers.csv`: `tests/parsers/bmg_mars/testdata/BMG Labtech format.CSV`
- `smart-control-fluorescence.xlsx`: `tests/parsers/bmg_labtech_smart_control/testdata/fluorescence_endpoint_example2.xlsx`
- `magellan-od-96.xlsx`: `tests/parsers/tecan_magellan/testdata/SLA01P01_OD.xlsx`
- `magellan-pro-compact.xlsx`: `tests/parsers/tecan_magellan/testdata/Magellan_Pro_7_5_compact_metadata.xlsx`
- `skanit-absorbance.xlsx`: `tests/parsers/thermo_skanit/testdata/skanit_data.xlsx`
- `skanit-luminescence-results-only.xlsx`: `tests/parsers/thermo_skanit/testdata/luciferase_assay_example.xlsx`
- `skanit-multi-plate-serial-date.xlsx`: `tests/parsers/thermo_skanit/testdata/multi-plate_example01.xlsx`

The expected decodes in `../plate_a-cases.json` were checked against the ASM
JSON allotropy keeps next to each input (software version, instrument model,
serial number, plate well count, detector and excitation wavelengths,
measurement time and measurement count).

The allotropy license, reproduced as it requires:

```
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
```
