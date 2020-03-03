from benchmark import Benchmark
from invocation import Invocation
from execution import Execution
from utility import *
from shutil import copyfile
import sys, importlib
import tmptool
import re

#
# Tooling for PRISM
#

# configuration
prism_bin = './fix-syntax ./prism'
prism_mem_args = '--javamaxmem 11g'

# instance specific settings

specific_settings = {}

def get_specific_setting(benchmark: Benchmark):
    """ returns the instance specific settings that were configured"""
    id = benchmark.get_identifier()
    if id not in specific_settings:
        return ''
    return specific_settings[id][0]

loaded = False
def assert_loaded():
    if not loaded:
        copyfile("tool.py", os.path.join(sys.path[0], "tmptool.py"))
        importlib.reload(sys.modules["tmptool"])

def get_name():
    """ should return the name of the tool as listed on http://qcomp.org/competition/2020/"""
    return "PRISM"

def is_benchmark_supported(benchmark : Benchmark):
    """returns True if the provided benchmark is supported by the tool and if the given benchmark should appear on the generated benchmark list"""

    if benchmark.is_prism():
        # check for unsupported property types
        if benchmark.is_reward_bounded_probabilistic_reachability() or benchmark.is_reward_bounded_expected_reward():
            return False
        return True

def get_prism_invocation_model_prop_instance(benchmark : Benchmark):
    args = "{} {} --property {}".format(benchmark.get_prism_program_filename(), benchmark.get_prism_property_filename(), benchmark.get_property_name())
    if benchmark.get_open_parameter_def_string() != "":
        args += " -const {}".format(benchmark.get_open_parameter_def_string())
    return args

def get_invocations(benchmark : Benchmark):
    """
    Returns a list of invocations that invoke the tool for the given benchmark.
    It can be assumed that the current directory is the directory from which execute_invocations.py is executed.
    For QCOMP 2020, this should return a list of invocations for all tracks in which the tool can take part. For each track an invocation with default settings has to be provided and in addition, an optimized setting (e.g., the fastest engine and/or solution technique for this benchmark) can be specified. Only information about the model type, the property type and the state space size are allowed to be used to tweak the parameters.
   
    If this benchmark is not supported, an empty list has to be returned.
    """

    if not is_benchmark_supported(benchmark):
        return []

    # Gather options that are needed for this particular benchmark for any invocation of PRISM
    benchmark_instance = get_prism_invocation_model_prop_instance(benchmark);
    
    invocations = []

    basic_args = "{}".format(prism_mem_args);

    # default settings
    default_args = basic_args
    default_inv = Invocation()
    default_inv.identifier = "default"
    default_inv.track_id = "epsilon-correct"
    default_inv.add_command("{} {} {}".format(prism_bin, default_args, benchmark_instance))
    invocations.append(default_inv)

    # specific settings                     !!!!only information about model type, property type and state space size via benchmark.get_num_states_tweak() may be used for tweaking
    specific_inv = Invocation()
    specific_inv.identifier = "specific"
    specific_inv.track_id = "epsilon-correct"
    specific_args = get_specific_setting(benchmark)
    specific_inv.add_command("{} {} {} {}".format(prism_bin, basic_args, benchmark_instance, specific_args))
    invocations.append(specific_inv)

    #### TODO: add default and specific invocations for other track_ids 'correct', 'probably-epsilon-correct', 'often-epsilon-correct', 'often-epsilon-correct-10-min'
    ### remember that different tracks have different precisions

    return invocations

def grep_for_result(benchmark : Benchmark, log) :
    # match 'Result: x (...)' and 'Result (...): x (...)' style results
    m = re.search('^Result.*: ([^( \n]+)', log, re.MULTILINE)
    if m is None:
        return m
    return m.group(1)

def get_result(benchmark : Benchmark, execution : Execution):
    """
    Returns the result of the given execution on the given benchmark.
    This method is called after executing the commands of the associated invocation.
    One can either find the result in the tooloutput (as done here) or
    read the result from a file that the tool has produced.
    The returned value should be either 'true', 'false', a decimal number, or a fraction.
    """
    invocation = execution.invocation
    log = execution.concatenate_logs()
    return grep_for_result(benchmark, log)
