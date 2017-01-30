#!/bin/bash

# basic model

../../bin/prism scheduler.prism scheduler.props -prop 2 -const prob1=0.25 -gridresolution 2 >| logs/scheduler_basic_time.0.25.gridres2.log

../../bin/prism scheduler.prism scheduler.props -prop 2 -const prob1=0.5 -gridresolution 2 >| logs/scheduler_basic_time.0.5.gridres2.log

../../bin/prism scheduler.prism scheduler.props -prop 2 -const prob1=0.75 -gridresolution 2 >| logs/scheduler_basic_time.0.75.gridres2.log

../../bin/prism scheduler.prism scheduler.props -prop 3 -const prob1=0.25 -gridresolution 2 >| logs/scheduler_basic_energy.0.25.gridres2.log

../../bin/prism scheduler.prism scheduler.props -prop 3 -const prob1=0.5 -gridresolution 2 >| logs/scheduler_basic_energy.0.5.gridres2.log

../../bin/prism scheduler.prism scheduler.props -prop 3 -const prob1=0.75 -gridresolution 2 >| logs/scheduler_basic_energy.0.75.gridres2.log

# random delay model

../../bin/prism scheduler_prob.prism scheduler.props -prop 2 -const prob1=0.25 -gridresolution 2 >| logs/scheduler_prob_time.0.25.gridres2.log

../../bin/prism scheduler_prob.prism scheduler.props -prop 2 -const prob1=0.5 -gridresolution 2 >| logs/scheduler_prob_time.0.5.gridres2.log

../../bin/prism scheduler_prob.prism scheduler.props -prop 2 -const prob1=0.75 -gridresolution 2 >| logs/scheduler_prob_time.0.75.gridres2.log

../../bin/prism scheduler_prob.prism scheduler.props -prop 3 -const prob1=0.25 -gridresolution 2 >| logs/scheduler_prob_energy.0.25.gridres2.log

../../bin/prism scheduler_prob.prism scheduler.props -prop 3 -const prob1=0.5 -gridresolution 2 >| logs/scheduler_prob_energy.0.5.gridres2.log

../../bin/prism scheduler_prob.prism scheduler.props -prop 3 -const prob1=0.75 -gridresolution 2 >| logs/scheduler_prob_energy.0.75.gridres2.log
