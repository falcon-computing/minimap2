# Minimap Flow
Falcon's version of minimap2

### Build
Use cmake to build the project. The following options control FPGA build settings.
* `USE_FPGA`: enable FPGA when OpenCL runtime is found. Default: `ON`
* `USE_BLAZE`: use Blaze integration for FPGA support. Default: `ON`
* `USE_NATIVE`: use native OpenCL runtime for FPGA support. Turn off `USE_BLAZE`. Default: `OFF`
* `USE_CKERNEL`: use C kernels in FPGA mode. Turn off `USE_BLAZE` and `USE_NATIVE`. Default: `OFF`

### Run
The following command shows the basic way to run mmap-flow.
```
~$ ./minimap-flow --x=<ALIGN_MODE>  \
                  --t=<NUM_THREADS>  \
                  --output_dir=<OUTPUT_DIR>  \
                  <REFERENCE_FILE>  \
                  [<READS_FILE1>] \
                  [<READS_FILE2>]
```
- `ALIGN_MODE`: value for `--x` option, indicating the alignment mode. Currently only support `"sr"`
- `NUM_THREADS`: number of threads for execution. Default: number of availiable CPUs on the platform
- `OUTPUT_DIR`: path to output directory for SAM/BAM output. Default: `""`
- `REFERENCE_FILE`: usually the minimizer reference file. If the original `".fasta"` file is provided, it will be converted into the new format (and dumped when given a dump file path).
- `READS_FILE*`: the reads file. If no reads file is provided, mmap-flow will only try to convert and dump the reference.

The following options are used with FPGA build and run.
- `--use_fpga`: enable FPGA execution in the FPGA build. Default: `yes`
- `--fpga_only`: disable CPU execution for FPGA-availabe stages. Default: `no`
- `--fpga_threads`: number for fpga workers. Default: `1`
- `--fpga_path`: path to FPGA bitstream. Only used in build with `USE_NATIVE`. Default: `""`
- `--blaze_conf`: path to Blaze MMSW task config file. Only available in build with `USE_BLAZE`. Default: `""`

These options are only available in the FPGA build.