//==============================================================================
//	
//	Copyright (c) 2014-
//	Authors:
//	* Xueyi Zou <xz972@york.ac.uk> (University of York)
//	
//------------------------------------------------------------------------------
//	
//	This file is part of PRISM.
//	
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//	
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//	
//==============================================================================

package explicit;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

import prism.PrismComponent;
import prism.PrismException;
import prism.PrismUtils;
import cern.colt.Arrays;
import explicit.rewards.MDPRewards;
import explicit.rewards.MDPRewardsSimple;

/**
 * Explicit-state model checker for Partially Observable Markov decision processes (POMDPs).
 */
public class POMDPModelChecker extends MDPModelChecker
{
	/**
	 * Create a new POMDPModelChecker, inherit basic state from parent (unless null).
	 */
	public POMDPModelChecker(PrismComponent parent) throws PrismException
	{
		super(parent);
	}

	/**
	 * Compute reachability/until probabilities.
	 * i.e. compute the min/max probability of reaching a state in {@code target},
	 * while remaining in those in @{code remain}.
	 * @param pomdp The POMDP
	 * @param remain Remain in these states (optional: null means "all")
	 * @param target Target states
	 * @param min Min or max probabilities (true=min, false=max)
	 */
	public ModelCheckerResult computeReachProbs(POMDP pomdp, BitSet target, boolean min) throws PrismException
	{
		ModelCheckerResult res = null;
		long timer;
		String stratFilename = null;

		// Local copy of setting
		POMDPSolnMethod pomdpSolnMethod = this.pomdpSolnMethod;

		// Switch to a supported method, if necessary
		if (!(pomdpSolnMethod == POMDPSolnMethod.RTDP_BEL || pomdpSolnMethod == POMDPSolnMethod.FIXED_GRID)) {
			pomdpSolnMethod = POMDPSolnMethod.FIXED_GRID;
			mainLog.printWarning("Switching to POMDP solution method \"" + pomdpSolnMethod.fullName() + "\"");
		}

		// Start probabilistic reachability
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting probabilistic reachability (" + (min ? "min" : "max") + ")...");

		// If required, create/initialise strategy storage
		if (genStrat || exportAdv) {
			stratFilename = exportAdvFilename;//"policyGraph.txt";
		}

		// Compute rewards
		switch (pomdpSolnMethod) {
		case RTDP_BEL:
			res = computeReachProbsRTDPBel(pomdp, target, min, stratFilename);
			break;
		case FIXED_GRID:
			res = computeReachProbsFixedGrid(pomdp, target, min, stratFilename);
			break;
		default:
			throw new PrismException("Unknown POMDP solution method " + pomdpSolnMethod.fullName());
		}

		// Finished probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.println("Probabilistic reachability took " + timer / 1000.0 + " seconds.");

		// Update time taken
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Compute expected reachability rewards using RTDP_Bel.
	 * Optionally, store optimal (memoryless) strategy info. 
	 * @param pomdp The POMDP
	 * @param mdpRewards The rewards
	 * @param target Target states
	 * @param min Min or max rewards (true=min, false=max)
	 * @param strat Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	protected ModelCheckerResult computeReachProbsRTDPBel(POMDP pomdp, BitSet target, boolean min, String stratFilename) throws PrismException
	{
		if (!(pomdp instanceof POMDPSimple)) {
			throw new PrismException("Sorry, RDTP_Bel does not support POMDP other than POMDPSimple.");
		}
		POMDPSimple simplePOMDP = (POMDPSimple) pomdp;

		ModelCheckerResult res;
		int iters;
		final int D = 10;
		boolean done = false;
		long timer;

		// Start real time dynamic programming (belief)
		timer = System.currentTimeMillis();
		mainLog.println("Starting real time dynamic programming (belief) (" + (min ? "min" : "max") + ")...");

		// Store num states
		int n = simplePOMDP.getNumStates();

		// Create solution vector(s)
		double soln[] = new double[n];

		//HashMap for storing real time values for discretized belief states, entries are allocated on the fly, values are updated every trial
		HashMap<DiscretizedBelief, Double> vhash = new HashMap<>();
		HashMap<DiscretizedBelief, Double> vhash_backUp = new HashMap<>();

		double[] initialBelief = simplePOMDP.getInitialBeliefInDist();
		mainLog.println("length of initialBelief: " + initialBelief.length);

		vhash.put(new DiscretizedBelief(initialBelief, D), 0.0);
		vhash_backUp.put(new DiscretizedBelief(initialBelief, D), 0.0);

		double[] targetBelief = new double[n];
		targetBelief[target.nextSetBit(0)] = 1.0;
		vhash.put(new DiscretizedBelief(targetBelief, D), 1.0);
		vhash_backUp.put(new DiscretizedBelief(targetBelief, D), 1.0);

		int numObservations = simplePOMDP.getNumObservations();
		mainLog.println("number of observations: " + numObservations);

		int numChoices;
		int chosenActionIndex; // the bold 'a'	
		double value, chosenValue;//chosenValue: Q(chosenActionIndex, b)

		// Start iterations
		iters = 0;
		//		maxIters=10;
		while (!done && iters < maxIters) {
			mainLog.println("\nValue Iteration: " + iters);
			double[] b = initialBelief;
			//sample state s with probability b(s);		
			int s = PrismUtils.pickFromDistribution(b, Math.random());
			while (!PrismUtils.isTargetBlief(b, target))//trial
			{
				mainLog.println("belief b: " + Arrays.toString(b));
				mainLog.println("state s: " + s + " values:" + simplePOMDP.statesList.get(s));

				//evaluate each action in b
				numChoices = simplePOMDP.getNumChoices(s);
				chosenValue = min ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
				chosenActionIndex = -1;
				for (int a = 0; a < numChoices; a++) {
					value = 0; // c(a,b)					
					for (int o = 0; o < numObservations; o++) {
						double observationProb = simplePOMDP.getObservationProbAfterAction(b, a, o);
						if (observationProb <= 1e-6) {
							continue;//this observation is impossible, go on with the next
						}
						double[] nextBelief = simplePOMDP.getBeliefInDistAfterActionAndObservation(b, a, o);
						DiscretizedBelief disretizedNextBelief = new DiscretizedBelief(nextBelief, D);
						if (vhash.get(disretizedNextBelief) == null) {
							vhash.put(disretizedNextBelief, 0.0);
						}
						value += observationProb * vhash.get(disretizedNextBelief);
					}
					mainLog.println("value for action " + a + ": " + value);

					//selction action that minimizes/maximizes Q(a,b), i.e. value
					if ((min && chosenValue - value > 1.0e-6) || (!min && value - chosenValue > 1.0e-6))//value<bestValue
					{
						chosenValue = value;
						chosenActionIndex = a;
					}
					if (Math.abs(value - chosenValue) < 1.0e-6)//value==chosenValue
					{
						//random tie broker
						chosenActionIndex = Math.random() < 0.5 ? a : chosenActionIndex;
					}
				}
				mainLog.println("chosen value: " + chosenValue);
				mainLog.println("chosen action index: " + chosenActionIndex);
				//update V(b) to the minimum Q(a,b), i.e. chosenValue
				vhash.put(new DiscretizedBelief(b, D), chosenValue);

				/********************sample next state with probability P_chosenActionIndex(nextState|s)*********************/
				Iterator<Entry<Integer, Double>> tranIter = simplePOMDP.getTransitionsIterator(s, chosenActionIndex);
				// cumalative distribution function
				ArrayList<Entry<Integer, Double>> cdf = new ArrayList<>();
				double sum = 0;
				while (tranIter.hasNext()) {
					Entry<Integer, Double> dist = tranIter.next();
					sum += dist.getValue();
					cdf.add(new AbstractMap.SimpleEntry<Integer, Double>(dist.getKey(), sum));
				}
				double prob = Math.random();
				int nextState = PrismUtils.pickFromDistribution(cdf, prob);
				mainLog.println("nextState s': " + nextState + " values:" + simplePOMDP.statesList.get(nextState));

				/********************sample observation o with probability Q_chosenActionIndex(o|nextState)*********************/
				int o = simplePOMDP.getObservation(nextState);
				mainLog.println("observation o: " + o + " values:" + simplePOMDP.observationsList.get(o));

				/********************compute next belief after action chosenActionIndex and observation o*********************/
				double[] nextBelief = simplePOMDP.getBeliefInDistAfterActionAndObservation(b, chosenActionIndex, o);
				mainLog.println("nextBelief: " + Arrays.toString(nextBelief));

				//go on with this trial
				mainLog.println("go on with this trial------------------------------------------------------>>> ");

				b = nextBelief;
				s = nextState;
			}
			// Check termination
			done = PrismUtils.hashMapsAreClose(vhash, vhash_backUp, termCritParam, termCrit == TermCrit.RELATIVE);
			// back up	
			Set<Map.Entry<DiscretizedBelief, Double>> entries = vhash.entrySet();
			mainLog.println("vhash:");
			for (Map.Entry<DiscretizedBelief, Double> entry : entries) {
				mainLog.println(entry.getKey() + "      " + entry.getValue());
				vhash_backUp.put(entry.getKey(), entry.getValue());
			}
			iters++;
			mainLog.println();
		}

		// Non-convergence is an error (usually)
		if (!done && errorOnNonConverge) {
			String msg = "Iterative method did not converge within " + iters + " iterations.";
			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
			throw new PrismException(msg);
		}

		double rewardsForInitialBelief = vhash.get(new DiscretizedBelief(initialBelief, D));
		for (int initialState : simplePOMDP.getInitialStates()) {
			soln[initialState] = rewardsForInitialBelief;

		}

		if (stratFilename != null) {
			throw new PrismException("Sorry, RDTP_Bel does not support adversory extraction.");
		}
		// Finished real tiem dynamic programming
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Real time dynamic programming (belief) (" + (min ? "min" : "max") + ")");
		mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Compute expected reachability rewards using Lovejoy's fixed-resolution grid approach.
	 * Optionally, store optimal (memoryless) strategy info. 
	 * @param pomdp The POMMDP
	 * @param mdpRewards The rewards
	 * @param target Target states
	 * @param inf States for which reward is infinite
	 * @param min Min or max rewards (true=min, false=max)
	 * @param strat Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	protected ModelCheckerResult computeReachProbsFixedGrid(POMDP pomdp, BitSet target, boolean min, String stratFilename) throws PrismException
	{
		if (!(pomdp instanceof POMDPSimple)) {
			throw new PrismException("Sorry, FixedGrid does not support POMDP other than POMDPSimple.");
		}
		POMDPSimple simplePOMDP = (POMDPSimple) pomdp;

		int theGridResolution = gridResolution;
		ModelCheckerResult res;
		int i, j, iters;
		boolean done = false;
		long timer;

		// Start fixed-resolution grid approximation
		timer = System.currentTimeMillis();
		mainLog.println("Starting fixed-resolution grid approximation (" + (min ? "min" : "max") + ")...");

		// Store num states
		int numStates = simplePOMDP.getNumStates();

		// Create solution vector(s)
		double soln[] = new double[numStates];

		//HashMap for storing real time values for the discretized grid belief states
		HashMap<Belief, Double> vhash = new HashMap<>();
		HashMap<Belief, Double> vhash_backUp = new HashMap<>();

		int numObservations = simplePOMDP.getNumObservations();
		int numUnobservations = simplePOMDP.getNumUnobservations();

		//find out the observations for the target states
		TreeSet<Integer> targetObservsSet = new TreeSet<>();
		for (int bit = target.nextSetBit(0); bit >= 0; bit = target.nextSetBit(bit + 1)) {
			targetObservsSet.add(simplePOMDP.getObservation(bit));

		}
		LinkedList<Integer> targetObservs = new LinkedList<>(targetObservsSet);

		//initialize the grid points
		ArrayList<Belief> gridPoints = new ArrayList<>();//the set of grid points (discretized believes)
		ArrayList<Belief> unknownGridPoints = new ArrayList<>();//the set of unknown grid points (discretized believes)
		ArrayList<ArrayList<Double>> assignment;
		boolean isTargetObserv;
		for (int so = 0; so < numObservations; so++) {
			ArrayList<Integer> unobservsForObserv = new ArrayList<>();
			for (int s = 0; s < numStates; s++) {
				if (so == simplePOMDP.getObservation(s)) {
					unobservsForObserv.add(simplePOMDP.getUnobservation(s));
				}
			}
			assignment = fullAssignment(unobservsForObserv.size(), theGridResolution);

			isTargetObserv = targetObservs.isEmpty() ? false : ((Integer) targetObservs.peekFirst() == so);
			if (isTargetObserv) {
				targetObservs.removeFirst();
			}

			for (ArrayList<Double> inner : assignment) {
				double[] bu = new double[numUnobservations];
				int k = 0;
				for (int unobservForObserv : unobservsForObserv) {
					bu[unobservForObserv] = inner.get(k);
					k++;
				}

				Belief g = new Belief(so, bu);
				gridPoints.add(g);
				if (isTargetObserv) {
					vhash.put(g, 1.0);
					vhash_backUp.put(g, 1.0);
				} else {
					unknownGridPoints.add(g);
					vhash.put(g, 0.0);
					vhash_backUp.put(g, 0.0);
				}
			}
		}

		int K = gridPoints.size();
		int unK = unknownGridPoints.size();
		mainLog.print("Grid statistics: resolution=" + theGridResolution);
		mainLog.println(", points=" + K + ", unknown points=" + unK);

		int numChoices;
		int chosenActionIndex; // the bold 'a'	
		double value, chosenValue;//chosenValue: Q(chosenActionIndex, b)

		//mainLog.println(pomdp);
		
		iters = 0;
		done = false;
		List<List<HashMap<Integer, Double>>> observationProbs = new ArrayList<>();//memoization for reuse
		List<List<HashMap<Integer, Belief>>> nextBelieves = new ArrayList<>();//memoization for reuse

		// Construct grid belief "MDP"
		// Iterate over all (unknown) grid points
		mainLog.println("Building belief space approximation...");
		for (i = 0; i < unK; i++) {
			Belief b = unknownGridPoints.get(i);
			double[] beliefInDist = b.toDistributionOverStates(pomdp);
			//mainLog.println("Belief " + i + ": " + b);
			//mainLog.print("Belief dist:");
			//mainLog.println(beliefInDist);
			numChoices = simplePOMDP.getNumChoicesForObservation(b.so);
			List<HashMap<Integer, Double>> action_observation_probs = new ArrayList<>();// for memoization
			List<HashMap<Integer, Belief>> action_observation_Believes = new ArrayList<>();// for memoization
			for (int a = 0; a < numChoices; a++) {
				//mainLog.println(i+"/"+unK+", "+a+"/"+numChoices);
				HashMap<Integer, Double> observation_probs = new HashMap<>();// for memoization
				HashMap<Integer, Belief> observation_believes = new HashMap<>();// for memoization
				simplePOMDP.computeObservationProbsAfterAction(beliefInDist, a, observation_probs);
				for (Map.Entry<Integer, Double> entry : observation_probs.entrySet()) {
					int o = entry.getKey();
					//mainLog.println(i+"/"+unK+", "+a+"/"+numChoices+", "+o+"/"+numObservations);
					Belief nextBelief = simplePOMDP.getBeliefAfterActionAndObservation(b, a, o);
					//mainLog.print(i + "/" + unK + ", " + a + "/" + numChoices + ", " + o + "/" + numObservations);
					//mainLog.println(" - " + entry.getValue() + ":" + nextBelief);
					observation_believes.put(o, nextBelief);
				}
				action_observation_probs.add(observation_probs);
				action_observation_Believes.add(observation_believes);
			}
			observationProbs.add(action_observation_probs);
			nextBelieves.add(action_observation_Believes);
		}
		
		// Start iterations
		mainLog.println("Solving belief space approximation...");
		while (!done && iters < maxIters) {
			//mainLog.println("iteration: " + iters);
			// Iterate over all (unknown) grid points
			for (i = 0; i < unK; i++) {
				Belief b = unknownGridPoints.get(i);
				numChoices = simplePOMDP.getNumChoicesForObservation(b.so);

				chosenValue = min ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
				chosenActionIndex = -1;
				for (int a = 0; a < numChoices; a++) {
					//mainLog.println(i+"/"+unK+", "+a+"/"+numChoices);
					value = 0; // c(a,b)

					for (Map.Entry<Integer, Double> entry : observationProbs.get(i).get(a).entrySet()) {
						int o = entry.getKey();
						//mainLog.println(i+"/"+unK+", "+a+"/"+numChoices+", "+o+"/"+numObservations);
						double observationProb = entry.getValue();
						Belief nextBelief = nextBelieves.get(i).get(a).get(o);

						//find discretized grid points to approximate the nextBelief
						ArrayList<double[]> subSimplex = new ArrayList<>();
						double[] lambdas = new double[nextBelief.bu.length];
						getSubSimplexAndLambdas(nextBelief.bu, subSimplex, lambdas, theGridResolution);

						//calculate the approximate value for the nextBelief
						double sum = 0;
						for (j = 0; j < lambdas.length; j++) {
							if (lambdas[j] >= 1e-6) {
								sum += lambdas[j] * vhash_backUp.get(new Belief(o, subSimplex.get(j)));
							}

						}
						value += observationProb * sum;
					}

					//select action that minimizes/maximizes Q(a,b), i.e. value
					if ((min && chosenValue - value > 1.0e-6) || (!min && value - chosenValue > 1.0e-6))//value<bestValue
					{
						chosenValue = value;
						chosenActionIndex = a;
					}
					if (Math.abs(value - chosenValue) < 1.0e-6)//value==chosenValue
					{ //random tie broker
						chosenActionIndex = Math.random() < 0.5 ? a : chosenActionIndex;
					}
				}
				//update V(b) to the chosenValue
				vhash.put(b, chosenValue);
			}

			// Check termination
			done = PrismUtils.hashMapsAreClose(vhash, vhash_backUp, termCritParam, termCrit == TermCrit.RELATIVE);
			// back up	
			Set<Map.Entry<Belief, Double>> entries = vhash.entrySet();
			for (Map.Entry<Belief, Double> entry : entries) {
				vhash_backUp.put(entry.getKey(), entry.getValue());
			}

			iters++;
		}
		// Non-convergence is an error (usually)
		if (!done && errorOnNonConverge) {
			String msg = "Iterative method did not converge within " + iters + " iterations.";
			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
			throw new PrismException(msg);
		}

		Belief initialBelief = simplePOMDP.getInitialBelief();
		//mainLog.println("initialBelief: " + initialBelief);
		//find discretized grid points to approximate the initialBelief
		ArrayList<double[]> subSimplex = new ArrayList<>();
		double[] lambdas = new double[initialBelief.bu.length];
		getSubSimplexAndLambdas(initialBelief.bu, subSimplex, lambdas, theGridResolution);

		//calculate the approximate value for the initialBelief
		double rewardsForInitialBelief = 0;
		for (j = 0; j < numUnobservations; j++) {
			if (lambdas[j] >= 1e-6) {
				rewardsForInitialBelief += lambdas[j] * vhash_backUp.get(new Belief(initialBelief.so, subSimplex.get(j)));
			}

		}

		for (int initialState : simplePOMDP.getInitialStates()) {
			soln[initialState] = rewardsForInitialBelief;
		}

		// extract optimal policy and store it in file named @code stratFilename
		MDPSimple mdp = new MDPSimple();
		BitSet mdpTarget = new BitSet();
		IndexedSet<Belief> exploredBelieves = new IndexedSet<>(true);
		LinkedList<Belief> toBeExploredBelives = new LinkedList<>();
		exploredBelieves.add(initialBelief);
		toBeExploredBelives.offer(initialBelief);
		mdp.addState();
		mdp.addInitialState(0);
		FileWriter policyGraphFileWriter = null;
		try {
			if (stratFilename != null) {
				File outFile = new File(stratFilename);
				policyGraphFileWriter = new FileWriter(outFile, false);
			}

			int src = -1;
			while (!toBeExploredBelives.isEmpty()) {
				Belief b = toBeExploredBelives.pollFirst();
				src++;
				if (PrismUtils.isTargetBlief(b.toDistributionOverStates(simplePOMDP), target)) {
					mdpTarget.set(src);
				}
				ArrayList<String> policies = extractBestActions(src, b, vhash, simplePOMDP, null, min, exploredBelieves, toBeExploredBelives, target, mdp);
				if (stratFilename != null)
					for (String policy : policies) {
						policyGraphFileWriter.write(policy);
					}
			}
			mainLog.println("\nStrategy DTMC:");
			mainLog.print(mdp.infoStringTable());
			//mainLog.println(mdp);
			//mainLog.println(mdpTarget);
			MDPModelChecker mcMDP = new MDPModelChecker(this);
			ModelCheckerResult mcRes = mcMDP.computeReachProbs(mdp, mdpTarget, true);
			mainLog.println("Outer bound: " + rewardsForInitialBelief);
			mainLog.println("Inner bound: " + mcRes.soln[0]);

		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (policyGraphFileWriter != null) {
				try {
					policyGraphFileWriter.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

		// Finished fixed-resolution grid approximation
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Fixed-resolution grid approximation (" + (min ? "min" : "max") + ")");
		mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Compute expected reachability rewards.
	 * i.e. compute the min/max reward accumulated to reach a state in {@code target}.
	 * @param pomdp The POMDP
	 * @param mdpRewards The rewards
	 * @param target Target states
	 * @param min Min or max rewards (true=min, false=max)
	 */
	public ModelCheckerResult computeReachRewards(POMDP pomdp, MDPRewards mdpRewards, BitSet target, boolean min) throws PrismException
	{
		ModelCheckerResult res = null;
		long timer;
		String stratFilename = null;
		// Local copy of setting
		POMDPSolnMethod pomdpSolnMethod = this.pomdpSolnMethod;

		// Switch to a supported method, if necessary
		if (!(pomdpSolnMethod == POMDPSolnMethod.RTDP_BEL || pomdpSolnMethod == POMDPSolnMethod.FIXED_GRID)) {
			pomdpSolnMethod = POMDPSolnMethod.FIXED_GRID;
			mainLog.printWarning("Switching to POMDP solution method \"" + pomdpSolnMethod.fullName() + "\"");
		}

		// Start expected reachability
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting expected reachability (" + (min ? "min" : "max") + ")...");

		// If required, create/initialise strategy storage
		if (genStrat || exportAdv) {
			stratFilename = exportAdvFilename;
		}

		// Compute rewards
		switch (pomdpSolnMethod) {
		case RTDP_BEL:
			res = computeReachRewardsRTDPBel(pomdp, mdpRewards, target, min, stratFilename);
			break;
		case FIXED_GRID:
			res = computeReachRewardsFixedGrid(pomdp, mdpRewards, target, min, stratFilename);
			break;
		default:
			throw new PrismException("Unknown POMDP solution method " + pomdpSolnMethod.fullName());
		}

		// Finished expected reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.println("Expected reachability took " + timer / 1000.0 + " seconds.");

		// Update time taken
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Compute expected reachability rewards using RTDP_Bel.
	 * Optionally, store optimal (memoryless) strategy info. 
	 * @param pomdp The POMDP
	 * @param mdpRewards The rewards
	 * @param target Target states
	 * @param min Min or max rewards (true=min, false=max)
	 * @param strat Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	protected ModelCheckerResult computeReachRewardsRTDPBel(POMDP pomdp, MDPRewards mdpRewards, BitSet target, boolean min, String stratFilename)
			throws PrismException
	{
		if (!(pomdp instanceof POMDPSimple)) {
			throw new PrismException("Sorry, RDTP_Bel does not support POMDP other than POMDPSimple.");
		}
		POMDPSimple simplePOMDP = (POMDPSimple) pomdp;

		ModelCheckerResult res;
		int iters;
		final int D = 10;
		boolean done = false;
		long timer;

		// Start real time dynamic programming (belief)
		timer = System.currentTimeMillis();
		mainLog.println("Starting real time dynamic programming (belief) (" + (min ? "min" : "max") + ")...");

		// Store num states
		int n = simplePOMDP.getNumStates();

		// Create solution vector(s)
		double soln[] = new double[n];

		//HashMap for storing real time values for discretized belief states, entries are allocated on the fly, values are updated every trial
		HashMap<DiscretizedBelief, Double> vhash = new HashMap<>();
		HashMap<DiscretizedBelief, Double> vhash_backUp = new HashMap<>();

		double[] initialBelief = simplePOMDP.getInitialBeliefInDist();
		//mainLog.println("length of initialBelief: " + initialBelief.length);

		vhash.put(new DiscretizedBelief(initialBelief, D), 0.0);
		vhash_backUp.put(new DiscretizedBelief(initialBelief, D), 0.0);

		int numObservations = simplePOMDP.getNumObservations();
		mainLog.println("number of observations: " + numObservations);

		int numChoices;
		int chosenActionIndex; // the bold 'a'	
		double value, chosenValue;//chosenValue: Q(chosenActionIndex, b)

		// Start iterations
		iters = 0;
		while (!done && iters < maxIters) {
			mainLog.println("\nValue Iteration: " + iters);
			double[] b = initialBelief;
			//sample state s with probability b(s);		
			int s = PrismUtils.pickFromDistribution(b, Math.random());
			while (!PrismUtils.isTargetBlief(b, target))//trial
			{
				mainLog.println("belief b: " + Arrays.toString(b));
				mainLog.println("state s: " + s + " values:" + simplePOMDP.statesList.get(s));

				//evaluate each action in b
				numChoices = simplePOMDP.getNumChoices(s);
				chosenValue = min ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
				chosenActionIndex = -1;
				for (int a = 0; a < numChoices; a++) {
					value = simplePOMDP.getCostAfterAction(b, a, mdpRewards); // c(a,b)

					for (int o = 0; o < numObservations; o++) {
						double observationProb = simplePOMDP.getObservationProbAfterAction(b, a, o);
						if (observationProb <= 1e-6) {
							continue;//this observation is impossible, go on with the next
						}
						double[] nextBelief = simplePOMDP.getBeliefInDistAfterActionAndObservation(b, a, o);
						DiscretizedBelief disretizedNextBelief = new DiscretizedBelief(nextBelief, D);
						if (vhash.get(disretizedNextBelief) == null) {
							vhash.put(disretizedNextBelief, 0.0);
						}
						value += observationProb * vhash.get(disretizedNextBelief);
					}
					mainLog.println("value for action " + a + ": " + value);

					//selction action that minimizes/maximizes Q(a,b), i.e. value
					if ((min && chosenValue - value > 1.0e-6) || (!min && value - chosenValue > 1.0e-6))//value<bestValue
					{
						chosenValue = value;
						chosenActionIndex = a;
					}
					if (Math.abs(value - chosenValue) < 1.0e-6)//value==chosenValue
					{
						//random tie broker
						chosenActionIndex = Math.random() < 0.5 ? a : chosenActionIndex;
					}
				}
				mainLog.println("chosen value: " + chosenValue);
				mainLog.println("chosen action index: " + chosenActionIndex);
				//update V(b) to the minimum Q(a,b), i.e. chosenValue
				vhash.put(new DiscretizedBelief(b, D), chosenValue);

				/********************sample next state with probability P_chosenActionIndex(nextState|s)*********************/
				Iterator<Entry<Integer, Double>> tranIter = simplePOMDP.getTransitionsIterator(s, chosenActionIndex);
				// cumalative distribution function
				ArrayList<Entry<Integer, Double>> cdf = new ArrayList<>();
				double sum = 0;
				while (tranIter.hasNext()) {
					Entry<Integer, Double> dist = tranIter.next();
					sum += dist.getValue();
					cdf.add(new AbstractMap.SimpleEntry<Integer, Double>(dist.getKey(), sum));
				}
				double prob = Math.random();
				int nextState = PrismUtils.pickFromDistribution(cdf, prob);
				mainLog.println("nextState s': " + nextState + " values:" + simplePOMDP.statesList.get(nextState));

				/********************sample observation o with probability Q_chosenActionIndex(o|nextState)*********************/
				int o = simplePOMDP.getObservation(nextState);
				mainLog.println("observation o: " + o + " values:" + simplePOMDP.observationsList.get(o));

				/********************compute next belief after action chosenActionIndex and observation o*********************/
				double[] nextBelief = simplePOMDP.getBeliefInDistAfterActionAndObservation(b, chosenActionIndex, o);
				mainLog.println("nextBelief: " + Arrays.toString(nextBelief));

				//go on with this trial
				mainLog.println("go on with this trial------------------------------------------------------>>> ");

				b = nextBelief;
				s = nextState;
			}
			// Check termination
			done = PrismUtils.hashMapsAreClose(vhash, vhash_backUp, termCritParam, termCrit == TermCrit.RELATIVE);
			// back up	
			Set<Map.Entry<DiscretizedBelief, Double>> entries = vhash.entrySet();
			mainLog.println("vhash:");
			for (Map.Entry<DiscretizedBelief, Double> entry : entries) {
				mainLog.println(entry.getKey() + "      " + entry.getValue());
				vhash_backUp.put(entry.getKey(), entry.getValue());
			}
			iters++;
			mainLog.println();
		}

		// Non-convergence is an error (usually)
		if (!done && errorOnNonConverge) {
			String msg = "Iterative method did not converge within " + iters + " iterations.";
			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
			throw new PrismException(msg);
		}

		double rewardsForInitialBelief = vhash.get(new DiscretizedBelief(initialBelief, D));
		for (int initialState : simplePOMDP.getInitialStates()) {
			soln[initialState] = rewardsForInitialBelief;
		}

		if (stratFilename != null) {
			throw new PrismException("Sorry, RDTP_Bel does not support adversory extraction.");
		}

		// Finished real tiem dynamic programming
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Real time dynamic programming (belief) (" + (min ? "min" : "max") + ")");
		mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	/**
	 * Compute expected reachability rewards using Lovejoy's fixed-resolution grid approach.
	 * Optionally, store optimal (memoryless) strategy info. 
	 * @param pomdp The POMMDP
	 * @param mdpRewards The rewards
	 * @param target Target states
	 * @param inf States for which reward is infinite
	 * @param min Min or max rewards (true=min, false=max)
	 * @param strat Storage for (memoryless) strategy choice indices (ignored if null)
	 */
	protected ModelCheckerResult computeReachRewardsFixedGrid(POMDP pomdp, MDPRewards mdpRewards, BitSet target, boolean min, String stratFilename)
			throws PrismException
	{
		if (!(pomdp instanceof POMDPSimple)) {
			throw new PrismException("Sorry, FixedGrid does not support POMDP other than POMDPSimple.");
		}
		POMDPSimple simplePOMDP = (POMDPSimple) pomdp;

		int theGridResolution = gridResolution;
		ModelCheckerResult res;
		int i, j, iters;
		boolean done = false;
		long timer;

		// Start fixed-resolution grid approximation
		timer = System.currentTimeMillis();
		mainLog.println("Starting fixed-resolution grid approximation (" + (min ? "min" : "max") + ")...");

		// Store num states
		int numStates = simplePOMDP.getNumStates();

		// Create solution vector(s)
		double soln[] = new double[numStates];

		//HashMap for storing real time values for the discretized grid belief states
		HashMap<Belief, Double> vhash = new HashMap<>();
		HashMap<Belief, Double> vhash_backUp = new HashMap<>();

		int numObservations = simplePOMDP.getNumObservations();
		int numUnobservations = simplePOMDP.getNumUnobservations();

		// find out the observations for the target states
		TreeSet<Integer> targetObservsSet = new TreeSet<>();
		for (int bit = target.nextSetBit(0); bit >= 0; bit = target.nextSetBit(bit + 1)) {
			targetObservsSet.add(simplePOMDP.getObservation(bit));
		}

		LinkedList<Integer> targetObservs = new LinkedList<>(targetObservsSet);

		//initialize the grid points
		ArrayList<Belief> gridPoints = new ArrayList<>();//the set of grid points (discretized believes)
		ArrayList<Belief> unknownGridPoints = new ArrayList<>();//the set of unknown grid points (discretized believes)
		ArrayList<ArrayList<Double>> assignment;
		boolean isTargetObserv;
		for (int so = 0; so < numObservations; so++) {
			ArrayList<Integer> unobservsForObserv = new ArrayList<>();
			for (int s = 0; s < numStates; s++) {
				if (so == simplePOMDP.getObservation(s)) {
					unobservsForObserv.add(simplePOMDP.getUnobservation(s));
				}
			}
			assignment = fullAssignment(unobservsForObserv.size(), theGridResolution);

			isTargetObserv = targetObservs.isEmpty() ? false : ((Integer) targetObservs.peekFirst() == so);
			if (isTargetObserv) {
				targetObservs.removeFirst();
			}

			for (ArrayList<Double> inner : assignment) {
				double[] bu = new double[numUnobservations];
				int k = 0;
				for (int unobservForObserv : unobservsForObserv) {
					bu[unobservForObserv] = inner.get(k);
					k++;
				}

				Belief g = new Belief(so, bu);
				gridPoints.add(g);
				if (isTargetObserv) {
					//					 mainLog.println(g);	
					vhash.put(g, 0.0);
					vhash_backUp.put(g, 0.0);
				} else {
					unknownGridPoints.add(g);
					vhash.put(g, 0.0);
					vhash_backUp.put(g, 0.0);
				}
			}
		}

		int K = gridPoints.size();
		int unK = unknownGridPoints.size();
		mainLog.print("Grid statistics: resolution=" + theGridResolution);
		mainLog.println(", points=" + K + ", unknown points=" + unK);

		int numChoices;
		int chosenActionIndex; // the bold 'a'	
		double value, chosenValue;//chosenValue: Q(chosenActionIndex, b)

		List<List<Double>> rewards = new ArrayList<>(); // using memoization for reuse
		List<List<HashMap<Integer, Double>>> observationProbs = new ArrayList<>();// using memoization for reuse
		List<List<HashMap<Integer, Belief>>> nextBelieves = new ArrayList<>();// using memoization for reuse
		iters = 0;
		done = false;

		// Construct grid belief "MDP"
		// Iterate over all (unknown) grid points
		mainLog.println("Building belief space approximation...");
		for (i = 0; i < unK; i++) {
			Belief b = unknownGridPoints.get(i);
			double[] beliefInDist = b.toDistributionOverStates(pomdp);
			//mainLog.println("Belief " + i + ": " + b);
			//mainLog.print("Belief dist:");
			//mainLog.println(beliefInDist);
			numChoices = simplePOMDP.getNumChoicesForObservation(b.so);
			List<Double> action_reward = new ArrayList<>();// for memoization
			List<HashMap<Integer, Double>> action_observation_probs = new ArrayList<>();// for memoization
			List<HashMap<Integer, Belief>> action_observation_Believes = new ArrayList<>();// for memoization
			for (int a = 0; a < numChoices; a++) {
				//mainLog.println(i+"/"+unK+", "+a+"/"+numChoices);
				action_reward.add(simplePOMDP.getCostAfterAction(b, a, mdpRewards)); // c(a,b)
				HashMap<Integer, Double> observation_probs = new HashMap<>();// for memoization
				HashMap<Integer, Belief> observation_believes = new HashMap<>();// for memoization
				simplePOMDP.computeObservationProbsAfterAction(beliefInDist, a, observation_probs);
				for (Map.Entry<Integer, Double> entry : observation_probs.entrySet()) {
					int o = entry.getKey();
					//mainLog.println(i+"/"+unK+", "+a+"/"+numChoices+", "+o+"/"+numObservations);
					Belief nextBelief = simplePOMDP.getBeliefAfterActionAndObservation(b, a, o);
					//mainLog.print(i + "/" + unK + ", " + a + "/" + numChoices + ", " + o + "/" + numObservations);
					//mainLog.println(" - " + entry.getValue() + ":" + nextBelief);
					observation_believes.put(o, nextBelief);
				}
				action_observation_probs.add(observation_probs);
				action_observation_Believes.add(observation_believes);
			}
			rewards.add(action_reward);
			observationProbs.add(action_observation_probs);
			nextBelieves.add(action_observation_Believes);
		}
		
		// Start iterations
		mainLog.println("Solving belief space approximation...");
		while (!done && iters < maxIters) {
			//mainLog.println("iteration: " + iters);
			// Iterate over all (unknown) grid points
			for (i = 0; i < unK; i++) {
				Belief b = unknownGridPoints.get(i);
				numChoices = simplePOMDP.getNumChoicesForObservation(b.so);

				chosenValue = min ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
				chosenActionIndex = -1;
				for (int a = 0; a < numChoices; a++) {
					value = rewards.get(i).get(a);

					for (Map.Entry<Integer, Double> entry : observationProbs.get(i).get(a).entrySet()) {
						int o = entry.getKey();
						//mainLog.println(i+"/"+unK+", "+a+"/"+numChoices+", "+o+"/"+numObservations);
						double observationProb = entry.getValue();
						Belief nextBelief = nextBelieves.get(i).get(a).get(o);
						
						//find discretized grid points to approximate the nextBelief
						ArrayList<double[]> subSimplex = new ArrayList<>();
						double[] lambdas = new double[nextBelief.bu.length];
						getSubSimplexAndLambdas(nextBelief.bu, subSimplex, lambdas, theGridResolution);

						//calculate the approximate value for the nextBelief
						double sum = 0;
						for (j = 0; j < lambdas.length; j++) {
							if (lambdas[j] >= 1e-6) {
								sum += lambdas[j] * vhash_backUp.get(new Belief(o, subSimplex.get(j)));
							}
						}
						value += observationProb * sum;
					}

					//select action that minimizes/maximizes Q(a,b), i.e. value
					if ((min && chosenValue - value > 1.0e-6) || (!min && value - chosenValue > 1.0e-6))//value<bestValue
					{
						chosenValue = value;
						chosenActionIndex = a;
					}
					if (Math.abs(value - chosenValue) < 1.0e-6)//value==chosenValue
					{
						//random tie broker
						chosenActionIndex = Math.random() < 0.5 ? a : chosenActionIndex;
					}
				}

				//update V(b) to the chosenValue
				vhash.put(b, chosenValue);
			}

			// Check termination
			done = PrismUtils.hashMapsAreClose(vhash, vhash_backUp, termCritParam, termCrit == TermCrit.RELATIVE);
			// back up	
			Set<Map.Entry<Belief, Double>> entries = vhash.entrySet();
			for (Map.Entry<Belief, Double> entry : entries) {
				vhash_backUp.put(entry.getKey(), entry.getValue());
			}

			iters++;
		}
		// Non-convergence is an error (usually)
		if (!done && errorOnNonConverge) {
			String msg = "Iterative method did not converge within " + iters + " iterations.";
			msg += "\nConsider using a different numerical method or increasing the maximum number of iterations";
			throw new PrismException(msg);
		}

		Belief initialBelief = simplePOMDP.getInitialBelief();
		mainLog.println("initialBelief: " + initialBelief);
		//find discretized grid points to approximate the initialBelief
		ArrayList<double[]> subSimplex = new ArrayList<>();
		double[] lambdas = new double[initialBelief.bu.length];
		getSubSimplexAndLambdas(initialBelief.bu, subSimplex, lambdas, theGridResolution);

		//calculate the approximate value for the initialBelief
		double rewardsForInitialBelief = 0;
		for (j = 0; j < numUnobservations; j++) {
			if (lambdas[j] >= 1e-6) {
				rewardsForInitialBelief += lambdas[j] * vhash_backUp.get(new Belief(initialBelief.so, subSimplex.get(j)));
			}

		}
		for (int initialState : simplePOMDP.getInitialStates()) {
			soln[initialState] = rewardsForInitialBelief;
		}

		// extract optimal adversary and store it to file named stratFilename
		MDPSimple mdp = new MDPSimple();
		BitSet mdpTarget = new BitSet();
		IndexedSet<Belief> exploredBelieves = new IndexedSet<>(true);
		LinkedList<Belief> toBeExploredBelives = new LinkedList<>();
		exploredBelieves.add(initialBelief);
		toBeExploredBelives.offer(initialBelief);
		mdp.addState();
		mdp.addInitialState(0);
		FileWriter policyGraphFileWriter = null;
		try {
			if (stratFilename != null) {
				File outFile = new File(stratFilename);
				policyGraphFileWriter = new FileWriter(outFile, false);
			}

			int src = -1;
			while (!toBeExploredBelives.isEmpty()) {
				Belief b = toBeExploredBelives.pollFirst();
				src++;
				if (PrismUtils.isTargetBlief(b.toDistributionOverStates(simplePOMDP), target)) {
					mdpTarget.set(src);
				}
				ArrayList<String> policies = extractBestActions(src, b, vhash, simplePOMDP, mdpRewards, min, exploredBelieves, toBeExploredBelives, target, mdp);
				if (stratFilename != null)
					for (String policy : policies) {
						policyGraphFileWriter.write(policy);
					}
			}
			mainLog.println("\nStrategy DTMC:");
			mainLog.print(mdp.infoStringTable());
			//mainLog.println(mdp);
			//mdp.exportToDotFile("adv.dot", mdpTarget);
			//mainLog.println(mdpTarget);
			MDPRewardsSimple mdpRewardsNew = new MDPRewardsSimple(mdp.getNumStates());
			for (int ii = 0; ii < mdp.getNumStates(); ii++) {
				Double rew = ((Double) mdp.getAction(ii, 0));
				mdpRewardsNew.addToStateReward(ii, rew == null ? 0.0 : rew.doubleValue());
			}
			MDPModelChecker mcMDP = new MDPModelChecker(this);
			ModelCheckerResult mcRes = mcMDP.computeReachRewards(mdp, mdpRewardsNew, mdpTarget, true);
			mainLog.println("Outer bound: " + rewardsForInitialBelief);
			mainLog.println("Inner bound: " + mcRes.soln[0]);

		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (policyGraphFileWriter != null) {
				try {
					policyGraphFileWriter.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

		// Finished fixed-resolution grid approximation
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Fixed-resolution grid approximation (" + (min ? "min" : "max") + ")");
		mainLog.println(" took " + iters + " iterations and " + timer / 1000.0 + " seconds.");

		// Return results
		res = new ModelCheckerResult();
		res.soln = soln;
		res.numIters = iters;
		res.timeTaken = timer / 1000.0;
		return res;
	}

	private ArrayList<ArrayList<Integer>> assignGPrime(int startIndex, int min, int max, int length)
	{
		ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
		if (startIndex == length - 1) {
			for (int i = min; i <= max; i++) {
				ArrayList<Integer> innerList = new ArrayList<>();
				innerList.add(i);
				result.add(innerList);
			}
		} else {
			for (int i = min; i <= max; i++) {
				ArrayList<ArrayList<Integer>> nextResult = assignGPrime(startIndex + 1, 0, i, length);
				for (ArrayList<Integer> nextReulstInner : nextResult) {
					ArrayList<Integer> innerList = new ArrayList<>();
					innerList.add(i);
					for (Integer a : nextReulstInner) {
						innerList.add(a);
					}
					result.add(innerList);
				}
			}
		}

		return result;
	}

	private ArrayList<ArrayList<Double>> fullAssignment(int length, int resolution)
	{
		ArrayList<ArrayList<Integer>> GPrime = assignGPrime(0, resolution, resolution, length);
		ArrayList<ArrayList<Double>> result = new ArrayList<ArrayList<Double>>();
		for (ArrayList<Integer> GPrimeInner : GPrime) {
			ArrayList<Double> innerList = new ArrayList<>();
			int i;
			for (i = 0; i < length - 1; i++) {
				int temp = GPrimeInner.get(i) - GPrimeInner.get(i + 1);
				innerList.add((double) temp / resolution);
			}
			innerList.add((double) GPrimeInner.get(i) / resolution);
			result.add(innerList);
		}
		return result;
	}

	private int[] getSortedPermutation(double[] inputArray)
	{
		int n = inputArray.length;
		double[] inputCopy = new double[n];
		int[] permutation = new int[n];
		int iState = 0, iIteration = 0;
		int iNonZeroEntry = 0, iZeroEntry = n - 1;
		boolean bDone = false;

		for (iState = n - 1; iState >= 0; iState--) {
			if (inputArray[iState] == 0.0) {
				inputCopy[iZeroEntry] = 0.0;
				permutation[iZeroEntry] = iState;
				iZeroEntry--;
			}

		}

		for (iState = 0; iState < n; iState++) {
			if (inputArray[iState] != 0.0) {
				inputCopy[iNonZeroEntry] = inputArray[iState];
				permutation[iNonZeroEntry] = iState;
				iNonZeroEntry++;
			}
		}

		while (!bDone) {
			bDone = true;
			for (iState = 0; iState < iNonZeroEntry - iIteration - 1; iState++) {
				if (inputCopy[iState] < inputCopy[iState + 1]) {
					swap(inputCopy, iState, iState + 1);
					swap(permutation, iState, iState + 1);
					bDone = false;
				}
			}
			iIteration++;
		}

		return permutation;
	}

	private void swap(int[] aiArray, int i, int j)
	{
		int temp = aiArray[i];
		aiArray[i] = aiArray[j];
		aiArray[j] = temp;
	}

	private void swap(double[] aiArray, int i, int j)
	{
		double temp = aiArray[i];
		aiArray[i] = aiArray[j];
		aiArray[j] = temp;
	}

	private boolean getSubSimplexAndLambdas(double[] b, ArrayList<double[]> subSimplex, double[] lambdas, int resolution)
	{
		int n = b.length;
		int M = resolution;

		double[] X = new double[n];
		int[] V = new int[n];
		double[] D = new double[n];
		for (int i = 0; i < n; i++) {
			X[i] = 0;
			for (int j = i; j < n; j++) {
				X[i] += M * b[j];
			}
			X[i] = Math.round(X[i] * 1e6) / 1e6;
			V[i] = (int) Math.floor(X[i]);
			D[i] = X[i] - V[i];
		}

		int[] P = getSortedPermutation(D);
		//		mainLog.println("X: "+ Arrays.toString(X));
		//		mainLog.println("V: "+ Arrays.toString(V));
		//		mainLog.println("D: "+ Arrays.toString(D));
		//		mainLog.println("P: "+ Arrays.toString(P));

		ArrayList<int[]> Qs = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			int[] Q = new int[n];
			if (i == 0) {
				for (int j = 0; j < n; j++) {
					Q[j] = V[j];
				}
				Qs.add(Q);
			} else {
				for (int j = 0; j < n; j++) {
					if (j == P[i - 1]) {
						Q[j] = Qs.get(i - 1)[j] + 1;
					} else {
						Q[j] = Qs.get(i - 1)[j];
					}

				}
				Qs.add(Q);
			}
			//			mainLog.println(Arrays.toString(Q));
		}

		for (int[] Q : Qs) {
			double[] node = new double[n];
			int i;
			for (i = 0; i < n - 1; i++) {
				int temp = Q[i] - Q[i + 1];
				node[i] = (double) temp / M;
			}
			node[i] = (double) Q[i] / M;
			subSimplex.add(node);
		}

		double sum = 0;
		for (int i = 1; i < n; i++) {
			double lambda = D[P[i - 1]] - D[P[i]];
			lambdas[i] = lambda;
			sum = sum + lambda;
		}
		lambdas[0] = 1 - sum;

		for (int i = 0; i < n; i++) {
			double sum2 = 0;
			for (int j = 0; j < n; j++) {
				sum2 += lambdas[j] * subSimplex.get(j)[i];
			}
			//			mainLog.println("b["+i+"]: "+b[i]+"  b^[i]:"+sum2);
			if (Math.abs(b[i] - sum2) > 1e-4) {
				return false;
			}

		}
		return true;
	}

	/**
	 * Find the best action for this belief state, add the belief state to the list
	 * of ones examined so far, and store the strategy info. We store this as an MDP.
	 * @param belief Belief state to examine
	 * @param vhash
	 * @param simplePOMDP
	 * @param mdpRewards
	 * @param min
	 * @param beliefList
	 */
	private ArrayList<String> extractBestActions(int src, Belief belief, HashMap<Belief, Double> vhash, POMDPSimple simplePOMDP, MDPRewards mdpRewards, boolean min,
			IndexedSet<Belief> exploredBelieves, LinkedList<Belief> toBeExploredBelives, BitSet target, MDPSimple mdp)
	{
		if (PrismUtils.isTargetBlief(belief.toDistributionOverStates(simplePOMDP), target)) {
			// Add self-loop
			/*Distribution distr = new Distribution();
			distr.set(src, 1);
			mdp.addActionLabelledChoice(src, distr, null);*/
			return new ArrayList<>();
		}
			
		ArrayList<String> policies = new ArrayList<>();
		int numObservations = simplePOMDP.getNumObservations();
		int RESOLUTION = gridResolution;

		double[] beliefInDist = belief.toDistributionOverStates(simplePOMDP);
		int numChoices = simplePOMDP.getNumChoicesForObservation(belief.so);
		double chosenValue = min ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
		int chosenActionIndex = -1;
		ArrayList<Integer> bestActions = new ArrayList<>();
		List<Double> action_reward = new ArrayList<>();
		List<HashMap<Integer, Double>> action_observation_probs = new ArrayList<>();
		List<HashMap<Integer, Belief>> action_observation_Believes = new ArrayList<>();
		//evaluate each action in b
		for (int a = 0; a < numChoices; a++) {
			double value = 0;
			if (mdpRewards != null) {
				value = simplePOMDP.getCostAfterAction(belief, a, mdpRewards); // c(a,b)	
			}
			// Build/store successor observations, probabilities and resulting beliefs
			HashMap<Integer, Double> observation_probs = new HashMap<>();
			HashMap<Integer, Belief> observation_believes = new HashMap<>();
			simplePOMDP.computeObservationProbsAfterAction(beliefInDist, a, observation_probs);
			for (Map.Entry<Integer, Double> entry : observation_probs.entrySet()) {
				int o = entry.getKey();
				Belief nextBelief = simplePOMDP.getBeliefAfterActionAndObservation(belief, a, o);
				observation_believes.put(o, nextBelief);
				double observationProb = observation_probs.get(o);

				//find discretized grid points to approximate the nextBelief
				ArrayList<double[]> subSimplex1 = new ArrayList<>();
				double[] lambdas1 = new double[nextBelief.bu.length];
				getSubSimplexAndLambdas(nextBelief.bu, subSimplex1, lambdas1, RESOLUTION);
				//calculate the approximate value for the nextBelief
				double sum = 0;
				for (int j = 0; j < lambdas1.length; j++) {
					if (lambdas1[j] >= 1e-6) {
						sum += lambdas1[j] * vhash.get(new Belief(o, subSimplex1.get(j)));
					}
				}
				value += observationProb * sum;
			}
			// Store the list of observations, probabilities and resulting beliefs for this action
			action_observation_probs.add(observation_probs);
			action_observation_Believes.add(observation_believes);

			//select action that minimizes/maximizes Q(a,b), i.e. value
			if ((min && chosenValue - value > 1.0e-6) || (!min && value - chosenValue > 1.0e-6))//value<bestValue
			{
				chosenValue = value;
				chosenActionIndex = a;
				bestActions.clear();
				bestActions.add(chosenActionIndex);
			} else if (Math.abs(value - chosenValue) < 1.0e-6)//value==chosenValue
			{
				//random tie broker
				chosenActionIndex = Math.random() < 0.5 ? a : chosenActionIndex;
				bestActions.clear();
				bestActions.add(a);
			}
		}

		Distribution distr = new Distribution();
		for (Integer a : bestActions) {
			for (Map.Entry<Integer, Double> entry : action_observation_probs.get(a).entrySet()) {
				int o = entry.getKey();
				double observationProb = entry.getValue();
				Belief nextBelief = action_observation_Believes.get(a).get(o);
				if (exploredBelieves.add(nextBelief)) {
					// If so, add to the explore list
					toBeExploredBelives.add(nextBelief);
					// And to model
					mdp.addState();
				}
				// Get index of state in state set
				int dest = exploredBelieves.getIndexOfLastAdd();

				String s = "";
				s += src + ":" + simplePOMDP.getObservationsList().get(belief.so) + "," + Arrays.toString(belief.bu);
				s += " & " + a + "/" + observationProb;
				s += " @ " + simplePOMDP.getObservationsList().get(o);
				s += " > " + dest + ":" + simplePOMDP.getObservationsList().get(nextBelief.so) + "," + Arrays.toString(nextBelief.bu);
				policies.add(s + "\n");
				//policies.add(Arrays.toString(belief.toDistributionOverStates(simplePOMDP))+"&"+a+"@"+simplePOMDP.getObservationsList().get(o)+">"+Arrays.toString(nextBelief.toDistributionOverStates(simplePOMDP))+"\n");
				distr.add(dest, observationProb);
			}
		}
		// Attach reward as action (hack)
		Object action = null;
		if (mdpRewards != null) {
			action = simplePOMDP.getCostAfterAction(belief, bestActions.get(0), mdpRewards);	
		}
		mdp.addActionLabelledChoice(src, distr, action);
		return policies;
	}

	/**
	 * Simple test program.
	 */
	public static void main(String args[])
	{
		POMDPModelChecker mc;
		POMDPSimple pomdp;
		ModelCheckerResult res;
		BitSet init, target;
		Map<String, BitSet> labels;
		boolean min = true;
		try {
			mc = new POMDPModelChecker(null);
			MDPSimple mdp = new MDPSimple();
			mdp.buildFromPrismExplicit(args[0]);
			//mainLog.println(mdp);
			labels = mc.loadLabelsFile(args[1]);
			//mainLog.println(labels);
			init = labels.get("init");
			target = labels.get(args[2]);
			if (target == null)
				throw new PrismException("Unknown label \"" + args[2] + "\"");
			for (int i = 3; i < args.length; i++) {
				if (args[i].equals("-min"))
					min = true;
				else if (args[i].equals("-max"))
					min = false;
				else if (args[i].equals("-nopre"))
					mc.setPrecomp(false);
			}
			pomdp = new POMDPSimple(mdp);
			res = mc.computeReachRewards(pomdp, null, target, min);
			System.out.println(res.soln[init.nextSetBit(0)]);
		} catch (PrismException e) {
			System.out.println(e);
		}
	}
}

class Tuple<X, Y, Z>
{
	public final X x;
	public final Y y;
	public final Z z;

	public Tuple(X x, Y y, Z z)
	{
		this.x = x;
		this.y = y;
		this.z = z;
	}

	@Override
	public String toString()
	{
		return (x + ", " + y + ", " + z);
	}
}