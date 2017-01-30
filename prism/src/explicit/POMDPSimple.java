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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import parser.Observation;
import parser.Unobservation;
import prism.ModelInfo;
import prism.ModelType;
import prism.PrismException;
import prism.PrismUtils;
import explicit.rewards.MDPRewards;

/**
 * Simple explicit-state representation of an POMDP.
 * The implementation is far from optimal, both in terms of memory usage and speed of access.
 * The model is, however, easy to manipulate. For a static model (i.e. one that does not change
 * after creation), consider MDPSparse, which is more efficient. 
 */
public class POMDPSimple extends MDPSimple implements POMDPExplicit
{
	/** Associated model info **/
	protected ModelInfo modelInfo;

	/** (Optionally) information about the observations of this model,
	 * i.e. the Observation object corresponding to each observation index. */
	public List<Observation> observationsList;

	/** (Optionally) information about the unobservations of this model,
	 * i.e. the Unobservation object corresponding to each unobservation index. */
	protected List<Unobservation> unobservationsList;

	// unobservable variables
	protected List<String> unobservableVars;

	/** Observable assigned to each state */
	protected List<Integer> observablesMap;

	// map a state to an unobservation
	public Map<Integer, Integer> state_unobserv_map;

	// Map from an observation to the number of possible choices
	public Map<Integer, Integer> observ_numChoices_map;

	protected int numObservations;
	protected int numUnobservations;

	// Constructors

	/**
	 * Constructor: empty POMDP.
	 */
	public POMDPSimple()
	{
		initialise(0);
	}

	/**
	 * Constructor: new POMDP with fixed number of states.
	 */
	public POMDPSimple(int numStates)
	{
		initialise(numStates);
	}

	/**
	 * Copy constructor.
	 */
	public POMDPSimple(POMDPSimple pomdp)
	{
		this(pomdp.numStates);
		copyFrom(pomdp);
		// Copy storage directly to avoid worrying about duplicate distributions (and for efficiency) 
		for (int s = 0; s < numStates; s++) {
			List<Distribution> distrs = trans.get(s);
			for (Distribution distr : pomdp.trans.get(s)) {
				distrs.add(new Distribution(distr));
			}
		}
		if (pomdp.actions != null) {
			actions = new ArrayList<List<Object>>(numStates);
			for (int s = 0; s < numStates; s++)
				actions.add(null);
			for (int s = 0; s < numStates; s++) {
				if (pomdp.actions.get(s) != null) {
					int n = pomdp.trans.get(s).size();
					List<Object> list = new ArrayList<Object>(n);
					for (int i = 0; i < n; i++) {
						list.add(pomdp.actions.get(s).get(i));
					}
					actions.set(s, list);
				}
			}
		}
		modelInfo = pomdp.modelInfo;
		observationsList = pomdp.observationsList;
		unobservationsList = pomdp.unobservationsList;
		observablesMap = pomdp.observablesMap;
		state_unobserv_map = pomdp.state_unobserv_map;
		observ_numChoices_map = pomdp.observ_numChoices_map;
		unobservableVars = pomdp.unobservableVars;
		numObservations = pomdp.numObservations;
		numUnobservations = pomdp.numUnobservations;

		// Copy flags/stats too
		allowDupes = pomdp.allowDupes;
		numDistrs = pomdp.numDistrs;
		numTransitions = pomdp.numTransitions;
		maxNumDistrs = pomdp.maxNumDistrs;
		maxNumDistrsOk = pomdp.maxNumDistrsOk;

		numObservations = pomdp.numObservations;
	}

	/**
	 * Construct an POMDP from an existing one and a state index permutation,
	 * i.e. in which state index i becomes index permut[i].
	 * Note: have to build new Distributions from scratch anyway to do this,
	 * so may as well provide this functionality as a constructor.
	 */
	public POMDPSimple(POMDPSimple pomdp, int permut[])
	{
		this(pomdp.numStates);
		copyFrom(pomdp, permut);
		// Copy storage directly to avoid worrying about duplicate distributions (and for efficiency)
		// (Since permut is a bijection, all structures and statistics are identical)
		for (int s = 0; s < numStates; s++) {
			List<Distribution> distrs = trans.get(permut[s]);
			for (Distribution distr : pomdp.trans.get(s)) {
				distrs.add(new Distribution(distr, permut));
			}
		}
		if (pomdp.actions != null) {
			actions = new ArrayList<List<Object>>(numStates);
			for (int s = 0; s < numStates; s++)
				actions.add(null);
			for (int s = 0; s < numStates; s++) {
				if (pomdp.actions.get(s) != null) {
					int n = pomdp.trans.get(s).size();
					List<Object> list = new ArrayList<Object>(n);
					for (int i = 0; i < n; i++) {
						list.add(pomdp.actions.get(s).get(i));
					}
					actions.set(permut[s], list);
				}
			}
		}
		modelInfo = pomdp.modelInfo;
		observationsList = pomdp.observationsList;
		unobservationsList = pomdp.unobservationsList;
		observablesMap = pomdp.observablesMap;
		state_unobserv_map = pomdp.state_unobserv_map;
		observ_numChoices_map = pomdp.observ_numChoices_map;
		unobservableVars = pomdp.unobservableVars;
		numObservations = pomdp.numObservations;
		numUnobservations = pomdp.numUnobservations;

		// Copy flags/stats too
		allowDupes = pomdp.allowDupes;
		numDistrs = pomdp.numDistrs;
		numTransitions = pomdp.numTransitions;
		maxNumDistrs = pomdp.maxNumDistrs;
		maxNumDistrsOk = pomdp.maxNumDistrsOk;
	}

	public POMDPSimple(MDPSimple mdp)
	{
		this(mdp.numStates);
		copyFrom(mdp);
		// Copy storage directly to avoid worrying about duplicate distributions (and for efficiency) 
		for (int s = 0; s < numStates; s++) {
			List<Distribution> distrs = trans.get(s);
			for (Distribution distr : mdp.trans.get(s)) {
				distrs.add(new Distribution(distr));
			}
			//observation of a state is the state itself
			observablesMap.set(s, s);
			//unobservation of a state null
			state_unobserv_map.put(s, null);
			//set the number of choices for an observation (state) 
			observ_numChoices_map.put(s, mdp.getNumChoices(s));
		}
		if (mdp.actions != null) {
			actions = new ArrayList<List<Object>>(numStates);
			for (int s = 0; s < numStates; s++)
				actions.add(null);
			for (int s = 0; s < numStates; s++) {
				if (mdp.actions.get(s) != null) {
					int n = mdp.trans.get(s).size();
					List<Object> list = new ArrayList<Object>(n);
					for (int i = 0; i < n; i++) {
						list.add(mdp.actions.get(s).get(i));
					}
					actions.set(s, list);
				}
			}
		}
		modelInfo = null;
		observationsList = null;
		unobservationsList = null;
		unobservableVars = null;
		numObservations = mdp.numStates;
		numUnobservations = 0;

		// Copy flags/stats too
		allowDupes = mdp.allowDupes;
		numDistrs = mdp.numDistrs;
		numTransitions = mdp.numTransitions;
		maxNumDistrs = mdp.maxNumDistrs;
		maxNumDistrsOk = mdp.maxNumDistrsOk;

		numObservations = mdp.numStates;
	}

	@Override
	public void addStates(int numToAdd)
	{
		super.addStates(numToAdd);
		if (observablesMap != null)
			for (int i = 0; i < numToAdd; i++) {
				observablesMap.add(-1);
		}
	}

	/**
	 * Set observation
	 */
	public void setObservation(int s, int o)
	{
		// If no observations array created yet, create it
		if (observablesMap == null) {
			observablesMap = new ArrayList<Integer>(numStates);
			for (int j = 0; j < numStates; j++)
				observablesMap.add(-1);
		}
		// Set observation
		observablesMap.set(s, o);
	}
	
	/**
	 * Set the associated model info.
	 */
	public void setModelInfo(ModelInfo modelInfo)
	{
		this.modelInfo = modelInfo;
	}

	/**
	 * Set the associated (read-only) observation list.
	 */
	public void setObservationsList(List<Observation> observationsList)
	{
		this.observationsList = observationsList;
	}

	/**
	 * Set the associated (read-only) observation list.
	 */
	public void setUnobservationsList(List<Unobservation> unobservationsList)
	{
		this.unobservationsList = unobservationsList;
	}

	/**
	 * Set the associated (read-only) unobservables list.
	 */
	public void setUnobservableVars(List<String> unobservableVars)
	{
		this.unobservableVars = unobservableVars;
	}

	/*************************** Mutators (for ModelSimple)********************************/

	@Override
	public void initialise(int numStates)
	{
		super.initialise(numStates);
		numDistrs = numTransitions = maxNumDistrs = 0;
		numObservations = 0;
		maxNumDistrsOk = true;
		trans = new ArrayList<List<Distribution>>(numStates);
		observablesMap = new ArrayList<Integer>();
		state_unobserv_map = new HashMap<Integer, Integer>();
		observ_numChoices_map = new HashMap<Integer, Integer>();
		for (int i = 0; i < numStates; i++) {
			trans.add(new ArrayList<Distribution>());
		}
		actions = null;
	}

	@Override
	public void clearState(int s)
	{
		// Do nothing if state does not exist
		if (s >= numStates || s < 0)
			return;
		// Clear data structures and update states
		List<Distribution> list = trans.get(s);
		numDistrs -= list.size();
		for (Distribution distr : list) {
			numTransitions -= distr.size();
		}
		maxNumDistrsOk = false;
		trans.get(s).clear();
		if (actions != null && actions.get(s) != null)
			actions.get(s).clear();
	}

	/*********************Accessors (for PartiallyObservableModel)*********************************************************/
	
	@Override
	public ModelInfo getModelInfo()
	{
		return modelInfo;
	}

	@Override
	public List<Observation> getObservationsList()
	{
		return observationsList;
	}

	@Override
	public List<Unobservation> getUnobservationsList()
	{
		return unobservationsList;
	}

	@Override
	public List<String> getUnobservableVars()
	{
		return unobservableVars;
	}

	/**
	 * Get the number of observations in state s.
	*/
	@Override
	public int getNumObservations()
	{
		return observationsList.size();
	}

	/**
	 * Get the number of unobservations in state s.
	*/
	@Override
	public int getNumUnobservations()
	{
		return unobservationsList.size();
	}

	/**
	 * Get the observation of state s.
	 */
	@Override
	public int getObservation(int s)
	{
		return observablesMap == null ? -1 : observablesMap.get(s);
	}

	/**
	 * Get the unobservation of state s.
	 */
	@Override
	public int getUnobservation(int s)
	{
		return state_unobserv_map.get(s);
	}

	/**
	* Get the probability of observing @code observation in state @code s.
	*/
	@Override
	public double getObservationProb(int s, int observation)
	{
		return this.getObservation(s) == observation ? 1.0 : 0.0;

	}

	/**
	 * Get the number of choices for an observation @code observ.
	 */
	@Override
	public int getNumChoicesForObservation(int observ)
	{
		return observ_numChoices_map.get(observ);
	}

	public void checkObservations() throws PrismException
	{
		// check if the same observation has the same choices
		observ_numChoices_map = new HashMap<Integer, Integer>();
		for (int s = 0; s < numStates; s++) {
			int numChoices = getNumChoices(s);
			Integer result = observ_numChoices_map.put(getObservation(s), numChoices);
			// mainLog.println("the result of adding it to map is "+result);
			if (result != null && result != numChoices) {
				throw new PrismException("For POMDP, the same observation should have the same actions. With the same observation "
						+ observationsList.get(getObservation(s)) + ", there may be " + result + " or " + numChoices + " choices!");
			}
		}
	}
	
	/*********************Accessors (for POMDP)*********************************************************/

	/**
	* Get initial belief state as an distribution over all states (array).
	*/
	@Override
	public Belief getInitialBelief()
	{
		double[] initialBeliefInDist = new double[numStates];
		for (Integer i : initialStates) {
			initialBeliefInDist[i] = 1;
		}
		PrismUtils.normalizeArray(initialBeliefInDist);
		return beliefInDistToBelief(initialBeliefInDist);
	}

	@Override
	public double[] getInitialBeliefInDist()
	{
		double[] initialBeliefInDist = new double[numStates];
		for (Integer i : initialStates) {
			initialBeliefInDist[i] = 1;
		}
		PrismUtils.normalizeArray(initialBeliefInDist);
		return initialBeliefInDist;
	}

	@Override
	public double getCostAfterAction(Belief belief, int action, MDPRewards mdpRewards)
	{
		double[] beliefInDist = belief.toDistributionOverStates(this);
		double cost = getCostAfterAction(beliefInDist, action, mdpRewards);
		return cost;
	}

	@Override
	public double getCostAfterAction(double[] beliefInDist, int action, MDPRewards mdpRewards)
	{
		double cost = 0;
		for (int i = 0; i < beliefInDist.length; i++) {
			if (beliefInDist[i] == 0) {
				cost += 0;
			} else {
				cost += beliefInDist[i] * (mdpRewards.getTransitionReward(i, action) + mdpRewards.getStateReward(i));
			}

		}
		return cost;
	}

	@Override
	public Belief getBeliefAfterAction(Belief belief, int action)
	{
		double[] beliefInDist = belief.toDistributionOverStates(this);
		double[] nextBeliefInDist = getBeliefInDistAfterAction(beliefInDist, action);
		return beliefInDistToBelief(nextBeliefInDist);
	}

	@Override
	public double[] getBeliefInDistAfterAction(double[] beliefInDist, int action)
	{
		int n = beliefInDist.length;
		double[] nextBeliefInDist = new double[n];
		for (int sp = 0; sp < n; sp++) {
			if (beliefInDist[sp] >= 1.0e-6) {
				Distribution distr = getChoice(sp, action);
				for (Map.Entry<Integer, Double> e : distr) {
					int s = (Integer) e.getKey();
					double prob = (Double) e.getValue();
					nextBeliefInDist[s] += beliefInDist[sp] * prob;
				}
			}
		}
		return nextBeliefInDist;
	}

	@Override // SLOW
	public double getObservationProbAfterAction(Belief belief, int action, int observation)
	{
		double[] beliefInDist = belief.toDistributionOverStates(this);
		double prob = getObservationProbAfterAction(beliefInDist, action, observation);
		return prob;
	}

	@Override // SLOW
	public double getObservationProbAfterAction(double[] beliefInDist, int action, int observation)
	{
		double[] beliefAfterAction = this.getBeliefInDistAfterAction(beliefInDist, action);
		int s;
		double prob = 0;
		for (s = 0; s < beliefAfterAction.length; s++) {
			prob += beliefAfterAction[s] * getObservationProb(s, observation);
		}
		return prob;
	}

	public void computeObservationProbsAfterAction(double[] beliefInDist, int action, HashMap<Integer, Double> observation_probs)
	{
		double[] beliefAfterAction = this.getBeliefInDistAfterAction(beliefInDist, action);
		for (int s = 0; s < beliefAfterAction.length; s++) {
			int o = getObservation(s);
			double probToAdd = beliefAfterAction[s];
			if (probToAdd > 1e-6) {
				Double lookup = observation_probs.get(o);
				if (lookup == null)
					observation_probs.put(o, probToAdd);
				else
					observation_probs.put(o, lookup + probToAdd);
			}
		}
	}
	
	@Override
	public Belief getBeliefAfterActionAndObservation(Belief belief, int action, int observation)
	{
		double[] beliefInDist = belief.toDistributionOverStates(this);
		double[] nextBeliefInDist = getBeliefInDistAfterActionAndObservation(beliefInDist, action, observation);
		Belief nextBelief = beliefInDistToBelief(nextBeliefInDist);
		if (nextBelief.so != observation) {
			System.err.println(nextBelief.so + "<--" + observation
					+ " something wrong with POMDPSimple.getBeliefAfterActionAndObservation(Belief belief, int action, int observation)");
		}
		return nextBelief;
	}

	@Override
	public double[] getBeliefInDistAfterActionAndObservation(double[] beliefInDist, int action, int observation)
	{
		int n = beliefInDist.length;
		double[] nextBelief = new double[n];
		double[] beliefAfterAction = this.getBeliefInDistAfterAction(beliefInDist, action);
		int i;
		double prob;
		for (i = 0; i < n; i++) {
			prob = beliefAfterAction[i] * getObservationProb(i, observation);
			nextBelief[i] = prob;
		}
		PrismUtils.normalizeArray(nextBelief);
		return nextBelief;
	}

	/************************ Standard methods*****************************************/

	@Override
	public String toString()
	{
		int i, j, n;
		Object o;
		String s = "";
		s = "[ ";
		for (i = 0; i < numStates; i++) {
			if (i > 0)
				s += ", ";
			s += i + "(" + getObservation(i) + "/" + getUnobservation(i) + "): ";
			s += "[";
			n = getNumChoices(i);
			for (j = 0; j < n; j++) {
				if (j > 0)
					s += ",";
				o = getAction(i, j);
				if (o != null)
					s += o + ":";
				s += trans.get(i).get(j);
			}
			s += "]";
		}
		s += " ]\n";
		return s;
	}

	@Override
	public boolean equals(Object o)
	{
		if (o == null || !(o instanceof POMDPSimple))
			return false;
		POMDPSimple mdp = (POMDPSimple) o;
		if (numStates != mdp.numStates)
			return false;
		if (!initialStates.equals(mdp.initialStates))
			return false;
		if (!trans.equals(mdp.trans))
			return false;
		// TODO: compare actions (complicated: null = null,null,null,...)
		return true;
	}

	@Override
	public ModelType getModelType()
	{
		return ModelType.POMDP;
	}

	@Override
	public String infoString()
	{
		String s = "";
		s += numStates + " states (" + getNumInitialStates() + " initial)";
		s += ", " + getNumTransitions() + " transitions";
		s += ", " + getNumChoices() + " choices";
		s += ", dist max/avg = " + getMaxNumChoices() + "/" + PrismUtils.formatDouble2dp(((double) getNumChoices()) / numStates);
		s += ", " + getNumObservations() + " observables";
		s += ", " + getNumUnobservations() + " unobservables";
		return s;
	}

	@Override
	public String infoStringTable()
	{
		String s = "";
		s += "States:        " + numStates + " (" + getNumInitialStates() + " initial)\n";
		s += "Transitions:   " + getNumTransitions() + "\n";
		s += "Choices:       " + getNumChoices() + "\n";
		s += "Max/avg:       " + getMaxNumChoices() + "/" + PrismUtils.formatDouble2dp(((double) getNumChoices()) / numStates) + "\n";
		s += "Observables:   " + getNumObservations() + "\n";
		s += "Unobservables: " + getNumUnobservations() + "\n";
		return s;
	}

	/***************************help methods*********************************************************/
	public Belief beliefInDistToBelief(double[] beliefInDist)
	{
		int so = -1;
		double[] bu = new double[this.getNumUnobservations()];
		for (int s = 0; s < beliefInDist.length; s++) {
			if (beliefInDist[s] != 0) {
				so = this.getObservation(s);
				bu[this.getUnobservation(s)] += beliefInDist[s];
			}
		}
		Belief belief = null;
		if (so != -1) {
			belief = new Belief(so, bu);
		} else {
			System.err.println("Something wrong in POMDPSimple.beliefInDistToBelief(double[] beliefInDist)");
		}
		return belief;
	}
}
