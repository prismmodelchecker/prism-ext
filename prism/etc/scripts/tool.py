from benchmark import Benchmark
from invocation import Invocation
from execution import Execution
import re

#
# Tooling for PRISM
#

# configuration
prism_bin = './fix-syntax ./prism'
prism_mem_args = '--javamaxmem 11g'

# instance specific settings

specific_settings = {
    'bluetooth.1.time': ('-mtbdd', ''),
    'crowds.5-20.positive': ('-hybrid', ''),
    'crowds.6-20.positive': ('-mtbdd',''),
    'egl.10-2.messagesB': ('-mtbdd',''),
    'egl.10-2.unfairA': ('-mtbdd',''),
    'egl.10-8.messagesB': ('-mtbdd',''),
    'egl.10-8.unfairA': ('-mtbdd',''),
    'haddad-monmege.100-0.7.exp_steps': ('-exact',''),
    'haddad-monmege.100-0.7.target': ('-exact',''),
    'herman.15.steps':  ('-sparse -gs','Select best engine and numerical method'),
    'nand.40-4.reliable': ('-sparse -bgs','Select best engine and backwards Gauss-Seidel as solution method, as the model is acyclic'),
    'nand.60-4.reliable': ('-sparse -bgs','Select best engine and backwards Gauss-Seidel as solution method, as the model is acyclic'),
    'oscillators.8-10-0.1-1-0.1-1.0.power_consumption': ('-hybrid',''),
    'oscillators.8-10-0.1-1-0.1-1.0.time_to_synch': ('-sparse -gs','Select best engine and numerical method'),
    'cluster.128-2000-20.premium_steady': ('-sparse -gs','Select best engine and numerical method'),
    'cluster.128-2000-20.qos1': ('-sparse',''),
    'cluster.64-2000-20.below_min': ('-sparse',''),
    'embedded.8-12.actuators': ('-hybrid',''),
    'embedded.8-12.up_time': ('-hybrid',''),
    'fms.8.productivity': ('-sparse -bgs','Select best engine and numerical method'),
    'kanban.5.throughput': ('-sparse -bgs', 'Select best engine and iteration method for CTMC steady-state'),
    'majority.2100.change_state': ('-sparse',''),
    'mapk_cascade.4-30.activated_time': ('-sparse -gs','Select best engine and numerical method'),
    'mapk_cascade.4-30.reactions': ('-sparse',''),
    'polling.18-16.s1_before_s2': ('-sparse -bgs','Select best engine and numerical method'),
    'speed-ind.2100.change_state': ('-sparse',''),
    'consensus.4-4.disagree': ('-sparse',''),
    'consensus.4-4.steps_min': ('-sparse',''),
    'consensus.6-2.disagree': ('-sparse -bgs','Select best engine and numerical method'),
    'consensus.6-2.steps_min': ('-sparse',''),
    'csma.3-4.all_before_max': ('-sparse',''),
    'csma.3-4.time_max': ('-hybrid',''),
    'csma.4-2.all_before_max': ('-hybrid',''),
    'csma.4-2.time_max': ('-hybrid',''),
#    'csma.3-4.some_before': ('',''),
    'eajs.5-250-11.ExpUtil': ('-sparse',''),
    'eajs.6-300-13.ExpUtil': ('-sparse',''),
    'pacman.100.crash': ('',''),
    'pacman.60.crash': ('',''),
    'pnueli-zuck.5.live': ('',''), # any symb engine is good, graph based
    'pnueli-zuck.10.live': ('',''), # any symb engine is good, graph based
    'rabin.10.live': ('',''),  # any symb engine is good, graph based
    'resource-gathering.1300-100-100.expgold': ('-sparse',''),
    'resource-gathering.1300-100-100.expsteps': ('-hybrid',''),
    'resource-gathering.1300-100-100.prgoldgem': ('-sparse',''),
    'wlan.4-0.sent': ('',''),  # any symb engine is good, graph based
    'wlan.4-0.cost_min': ('-sparse',''),
    'wlan.5-0.sent': ('',''),  # any symb engine is good, graph based
    'wlan.5-0.cost_min': ('-sparse',''),
    'wlan.6-0.sent': ('',''),  # any symb engine is good, graph based
    'wlan.6-0.cost_min': ('-sparse',''),
    'zeroconf.1000-8-false.correct_max': ('-explicit',''),
    'zeroconf.1000-8-false.correct_min': ('-explicit',''),
    'firewire-pta.30-5000.deadline': ('',''),  # PTA: use default STPG based analysis
    'firewire-pta.30-5000.eventually': ('',''),  # PTA: use default STPG based analysis
    'repudiation_malicious.20.deadline': ('',''),  # PTA: use default STPG based analysis
    'repudiation_malicious.20.eventually': ('',''),  # PTA: use default STPG based analysis
    'zeroconf-pta.200.deadline': ('',''),  # PTA: use default STPG based analysis
    'zeroconf-pta.200.incorrect': ('',''),  # PTA: use default STPG based analysis
}

def get_specific_setting(benchmark: Benchmark):
    """ returns the instance specific settings that were configured"""
    id = benchmark.get_identifier()
    if id not in specific_settings:
        return ''
    return specific_settings[id][0]

def get_specific_setting_explanation(benchmark: Benchmark):
    """ returns the explanation for the instance specific settings that were configured"""
    id = benchmark.get_identifier()
    if id not in specific_settings:
        return ''
    ex = specific_settings[id][1]
    if ex == '':
        if specific_settings[id][0] == '':
            ex = 'Use default settings'
        else:
            ex = 'Select best engine'
    return ex

def get_name():
    """ should return the name of the tool as listed on http://qcomp.org/competition/2019/"""
    return "PRISM"

def is_benchmark_supported(benchmark : Benchmark):
    """ Auxiliary function that returns True if the provided benchmark is supported by the tool"""

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
    For QCOMP 2019, this should return a list of size at most two, where
    the first entry (if present) corresponds to the default configuration of the tool and
    the second entry (if present) corresponds to an optimized setting (e.g., the fastest engine and/or solution technique for this benchmark).
    Please only provide two invocations if there is actually a difference between them.
    If this benchmark is not supported, an empty list has to be returned.
    For testing purposes, the script also allows to return more than two invocations.
    """

    if not is_benchmark_supported(benchmark):
        return []

    # Gather options that are needed for this particular benchmark for any invocation of Storm
    benchmark_instance = get_prism_invocation_model_prop_instance(benchmark);

    invocations = []


    basic_args = "{}".format(prism_mem_args);

    # default settings

    default_args = basic_args

    default_inv = Invocation()
    default_inv.identifier = "default"
    default_inv.note = "Default settings."
    default_inv.add_command("{} {} {}".format(prism_bin, default_args, benchmark_instance))
    invocations.append(default_inv)

    # specific settings
    specific_inv = Invocation()
    specific_inv.identifier = "specific"
    specific_args = get_specific_setting(benchmark)
    specific_inv.note = get_specific_setting_explanation(benchmark)
    specific_inv.add_command("{} {} {} {}".format(prism_bin, basic_args, benchmark_instance, specific_args))
    invocations.append(specific_inv)

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
