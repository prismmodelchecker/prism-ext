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

import java.util.List;

import parser.Observation;
import parser.Unobservation;
import prism.ModelInfo;

/**
 * Interface for (abstract) classes that provide (read-only) access to an explicit-state model with Partial Observability.
 */
public interface PartiallyObservableModel extends Model
{
	// Accessors
	/**
	 * Get the associated model info.
	 */
	public ModelInfo getModelInfo();

	/**
	 * Get access to a list of observations (optionally stored).
	 */
	public List<Observation> getObservationsList();

	/**
	 * Get access to a list of unobservations (optionally stored).
	 */
	public List<Unobservation> getUnobservationsList();

	/**
	 * Get access to a list of unobservable varialbes (optionally stored).
	 */
	public List<String> getUnobservableVars();

	/**
	 * Get the total number of observations over all states.
	 */
	public int getNumObservations();

	/**
	 * Get the total number of unobservations over all states.
	 */
	public int getNumUnobservations();

	/**
	 * Get the observation of state @codes.
	 */
	public int getObservation(int s);

	/**
	 * Get the unobservation of state @codes.
	 */
	public int getUnobservation(int s);

	/**
	 * Get the probability of observing @code observation in state @code s.
	 */
	public double getObservationProb(int s, int observation);

	/**
	 * Get the number of choices for an observation @code observ.
	 */
	public int getNumChoicesForObservation(int observ);
}