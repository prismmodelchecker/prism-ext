#!/bin/bash

# expected reachability

../../bin/prism pump.prism pump.props -prop 2 -const h0=2,h1=16,N=2 -gridresolution 2 >| logs/pump2.16.2.gridres2.log

../../bin/prism pump.prism pump.props -prop 2 -const h0=2,h1=16,N=2 -gridresolution 40 >| logs/pump2.16.2.gridres40.log

../../bin/prism pump.prism pump.props -prop 2 -const h0=2,h1=16,N=16 -gridresolution 2 >| logs/pump2.16.16.gridres2.log

# time-bounded probabilistic reachability

../../bin/prism pump_deadline.prism pump.props -prop 2 -const h0=2,h1=8,N=4,D=50 -gridresolution 2 >| logs/pump_deadline2.8.4.50.gridres2.log

../../bin/prism pump_deadline.prism pump.props -prop 2 -const h0=2,h1=8,N=4,D=50 -gridresolution 12 >| logs/pump_deadline2.8.4.50.gridres12.log

../../bin/prism pump_deadline.prism pump.props -prop 2 -const h0=2,h1=16,N=8,D=50 -gridresolution 2 >| logs/pump_deadline2.16.8.50.gridres2.log

../../bin/prism pump_deadline.prism pump.props -prop 2 -const h0=2,h1=16,N=8,D=100 -gridresolution 2 >| logs/pump_deadline2.16.8.100.gridres2.log
