## Wavefront DTW FPGA Accelerator (Jan 2026 – Feb 2026)

Designed and implemented a dedicated Wavefront Dynamic Time Warping (DTW) accelerator on the Xilinx Alveo U55C platform.

### Target Board
- Xilinx Alveo U55C

### Optimization Techniques
- Loop pipelining
- Loop unrolling
- Array partitioning
- Shifted buffer design for data reuse
- Bit-width optimization using ap_(u)int
- fully using bus width with vector

### Performance
- 15× speedup over sequential CPU implementation
- 136× speedup over naïve wavefront implementation
- Achieved Initiation Interval (II) = 1 in the core inner computation loop