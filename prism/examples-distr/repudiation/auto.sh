#!/bin/bash

# basic model

../../bin/prism repudiation.prism repudiation.props -const K=4 -gridresolution 8 >| logs/repudiation_basic4_gridres8.log

../../bin/prism repudiation.prism repudiation.props -const K=4 -gridresolution 24 >| logs/repudiation_basic4_gridres24.log

../../bin/prism repudiation.prism repudiation.props -const K=8 -gridresolution 4 >| logs/repudiation_basic8_gridres4.log

../../bin/prism repudiation.prism repudiation.props -const K=8 -gridresolution 8 >| logs/repudiation_basic8_gridres8.log

# complex model

../../bin/prism repudiation_complex.prism repudiation.props -const K=4 -gridresolution 4 >| logs/repudiation_complex4_gridres4.log

../../bin/prism repudiation_complex.prism repudiation.props -const K=4 -gridresolution 12 >| logs/repudiation_complex4_gridres12.log

../../bin/prism repudiation_complex.prism repudiation.props -const K=8 -gridresolution 2 >| logs/repudiation_complex8_gridres2.log

../../bin/prism repudiation_complex.prism repudiation.props -const K=8 -gridresolution 4 >| logs/repudiation_complex8_gridres4.log
