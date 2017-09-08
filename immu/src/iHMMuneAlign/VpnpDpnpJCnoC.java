package iHMMuneAlign;

/**
 * Rewrite by Jessica Ye
 * @author harris
 */

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.biojava.bio.Annotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dp.DotState;
import org.biojava.bio.dp.EmissionState;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.SimpleDotState;
import org.biojava.bio.dp.SimpleEmissionState;
import org.biojava.bio.dp.SimpleMarkovModel;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;


/**
 * this class funtions like a factory Instance, and will create a MarkovModel
 * for the Prediction of the CDR3 region in the Immuno Globulin System.
 * 
 * Every method other than "createModel(...)" is only a convenience method, a
 * collection of statements necessary to produce a Markov Model with the desired
 * functionality
 */
public class VpnpDpnpJCnoC implements GlobalDefines {


	// Gene type, V,D or J
	final int GENE_TYPE_V = 1;
    final int GENE_TYPE_D = 2;
	final int GENE_TYPE_J = 3;

	// globals
	boolean DEBUG = false; // flag for debugging

	// minimum exo nuclease probability before array is cut off
	// basically, the exo nuclease probability values are given
	// as a "mean" and a "standard deviation" value for each gene family
	// From these values, for each family, we create an array with the
	// exo. nucl. removal probabilites from zero removals until
	// the length of the entire gene or until the probability is
	// less than the "MIN_PROB_LIMIT" below, whichever comes first
	double MIN_PROB_LIMIT = 0.000001;

	// /////////////// Some Variables were made global to reduce number of
	// parameters passed in methods

	MutationSpectrum G_mutationSpectrum; // an instance of the mutation spectrum class

	// ////// two different ways of calculating the mutability scores
	MutabilityScore G_mutability_score; // container of methods for calculating
										// the mutability of a pentanucleotide

	ExponentialDecay G_exponential_decay; // container of methods for
										  // calculating the exponential decay
										  // in probability for a nucl. in a
										  // gene

	PrintWriter G_fostream; // output important data if needed

	/**
	 * constructor for this class
	 */
	public VpnpDpnpJCnoC(MutabilityScore mutabilityScore) {
		
		// assign the mutability scoring object to be used for the life span 
		G_mutability_score = mutabilityScore;

		// the object for calculating the exponential decay
	   G_exponential_decay = new ExponentialDecay();

		// make a MutationSpectrum lookup object from file
		this.G_mutationSpectrum = new MutationSpectrum("/srvr/ihmmune/tomcat7/webapps/immuFile/Mutation spectrum.txt");

		try {
			// this file stream is used for writing the nucleotide mutation
			// probabilites to a file. The output data format is a little cryptic
			G_fostream = new PrintWriter(new FileWriter(
					"/srvr/ihmmune/tomcat7/webapps/immuFile/gene nucleotide mutation probabilities.txt"));
		} catch (IOException ioex) {
			throw new Error("output stream could not be opened");
		}//--catch()

	}//--constructor

	/**
	 * this is the only public method of this class and will create the
	 * MarkovModel It starts the creation of a markov model with states
	 * Magical->V->N1->Dx->N2->Jx-N3->Magical Explanation of Variables
	 * "vSequence" : a VGene sequence, which is already pre aligned with the
	 * CDR3 sequence "VGene_start_offset": the number of nts at the beginning of
	 * the VGene which were not part of the pre alignment "completeVGeneLength":
	 * simply the length of the pre aligned VGene "probHolder": an instance of
	 * class ProbabilityHolder holding all the probabilities needed to create
	 * this MarkovModel "removed_C_Region": a flag indicating whether the CDR 3
	 * sequence had any C-Region (if it had, it is now removed) 'A_probability":
	 * the probability at position ZERO (VGene position ZERO) for the
	 * exponentially decaying mutation probability for a CDR 3 sequence
	 */
	public MarkovModel createModel(Sequence vSequence, int VGene_start_offset,
			int completeVGeneLength, File D_GENES_FILE, Sequence jSequence,
			ProbabilityHolder probHolder, boolean removed_C_Region,
			double A_probability) {
		//set up the dice alphabet
		FiniteAlphabet dna = DNATools.getDNA();
		FastaReader fr = new FastaReader();

		final int[] advance = { 1 }; // how many states to advance when going
									 // through model (always one for us)

		// List of V,D and J Gene states
		ArrayList vStates = new ArrayList();
		ArrayList dStates = new ArrayList(); // list of lists
		ArrayList jStates = new ArrayList(); // list of lists

		// list of Gene Sequences
		ArrayList dSequences = new ArrayList();
		ArrayList jSequences = new ArrayList();

		// list of lists of distributions for all genes
		ArrayList vStateDistributions = new ArrayList(); // not a list of lists, one list only
		ArrayList dStateDistributions = new ArrayList(); // list of lists
		ArrayList jStateDistributions = new ArrayList(); // list of lists

		// read the sequences from files in fasta format
		try {

			dSequences = fr.readFile(D_GENES_FILE);
			jSequences.add(jSequence);

		} catch (Exception e) {
			System.out.println("fr.readFile(): " + e.getMessage()
					+ D_GENES_FILE + " not a fasta file ???");
			System.exit(0);
		}// --catch

		// /////////////////// Creating the Markov Model /////////////////// //
		// In order to create this Markov Model without any problems, we must do
		// it in four separate
		// actions:
		// 1: create all the states in the model
		// 2: add all the created states to the model
		// 3: create all the state transitions in the model
		// 4: set all the state transition probabilities in the model

		//*************** create States *****************//

		// Dot States are the binding material between emission states when
		// we want to simplify and reduce the number of connections/transitions
		// between
		// states in the model (it looks better and it reduces computation time

		// create Dot States for use in this model
		// some connections require arrays of Dot State, because they are
		// connecting multiple
		// D / J genes to something else (say for example there is one X6 state
		// for each D gene)
		DotState X1a, X1b, X2, X3, X8, X9; // single dot states
		DotState[] X5b_states, X5a_states, X6_states, X7a_states, X7b_states, X11b_states, X11a_states, X12_states;

		// create the dot states and give them a name
		X1a = createDotState("X1a");
		X1b = createDotState("X1b");
		X2 = createDotState("X2");
		X3 = createDotState("X3");
		X8 = createDotState("X8");
		X9 = createDotState("X9");

		// create the dot states that are in fact multiple states in one
		// (arrays)
		X5a_states = createDotStateArray("X5a", dSequences.size());
		X5b_states = createDotStateArray("X5b", dSequences.size());
		X6_states = createDotStateArray("X6", dSequences.size());
		X7a_states = createDotStateArray("X7a", dSequences.size());
		X7b_states = createDotStateArray("X7a", dSequences.size());

		X11a_states = createDotStateArray("X11a", jSequences.size());
		X11b_states = createDotStateArray("X11b", jSequences.size());
		X12_states = createDotStateArray("X12", jSequences.size());

		// create V, D and J states
		createVStates(advance, dna, vSequence, vStates, vStateDistributions,
				VGene_start_offset, completeVGeneLength, probHolder,
				A_probability);
		createDStates(advance, dna, dSequences, dStates, dStateDistributions,
				probHolder, completeVGeneLength, A_probability);
		createJStates(advance, dna, jSequences, jStates, jStateDistributions,
				probHolder, completeVGeneLength, A_probability);

		// create the N states ( size must equal the information given in the probHolder data )
		EmissionState[] VD_N_states = createNStates(advance, dna,
				"V-D N-Region", true, probHolder.VD_N_Fields.length);
		
		EmissionState[] DJ_N_states = createNStates(advance, dna,
				"D-J N-Region", false, probHolder.DJ_N_Fields.length);
	

		// create the P - states
		// some of the P-regions (P2,P3,P4) are in fact multiple P-regions 
		// The exception of multiple P-regions for a P-region is the End of the
		// VGene P-region, because there is only one VGene

		// create single P-region
		EmissionState[] V_P_end_states = createPStateArray(advance, dna,
				"V-End P-Region", probHolder.V_end_P_Fields.length,
				V_END_P_TOKEN, vSequence, false);

		// create multiple P-regions
		ArrayList D_P_start_array_list = createPStatesArrayList(advance, dna,
				"D-Start P_Region", probHolder.D_start_P_Fields.length,
				D_START_P_TOKEN, dSequences, true);
		ArrayList D_P_end_array_list = createPStatesArrayList(advance, dna,
				"D-End P_Region", probHolder.D_end_P_Fields.length,
				D_END_P_TOKEN, dSequences, false);
		ArrayList J_P_start_array_list = createPStatesArrayList(advance, dna,
				"J-Start P_Region", probHolder.J_start_P_Fields.length,
				J_START_P_TOKEN, jSequences, true);

		//*************** create the MarkovModel *****************//

		// make an instance of a SimpleMarkovModel and refer to it as "vdj"
		// this is the model to which all the states are added
		SimpleMarkovModel vdj = new SimpleMarkovModel(1, dna, "VDJGene");

		//*************** add the states to the MarkovModel *****************//

		// add the Dot states to this markov model
		addDotStates(vdj, X1a, X1b, X2, X3, X8, X9);

		addDotStateArray(vdj, X5b_states);
		addDotStateArray(vdj, X5a_states);
		addDotStateArray(vdj, X6_states);

		addDotStateArray(vdj, X7a_states);
		addDotStateArray(vdj, X7b_states);

		addDotStateArray(vdj, X11b_states);
		addDotStateArray(vdj, X11a_states);
		addDotStateArray(vdj, X12_states);

		// add the P states
		addPStates(vdj, V_P_end_states);

		addPStates_list(vdj, D_P_start_array_list);
		addPStates_list(vdj, D_P_end_array_list);
		addPStates_list(vdj, J_P_start_array_list);

		// add the N states to this markov model
		addNStates(VD_N_states, vdj);
		addNStates(DJ_N_states, vdj);

		// add the V D J states to this markov model
		addVStates(vdj, vStates);
		addDStates(vdj, dStates);
		addJStates(vdj, jStates);

		//*************** create the state transitions *****************//

		// all state transitions must be created and added to the model,
		// then we can set the transition probabilities
		createVTransitions(vdj, vStates, X1a, X1b, probHolder);

		// create P transitions
		create_P_transitions(vdj, V_P_end_states, X2);

		create_P_transitions_multi_start(vdj, X5a_states, D_P_start_array_list,
				X6_states);
		create_P_transitions_multi_end(vdj, D_P_end_array_list, X8);
		create_P_transitions_multi_start(vdj, X11a_states,
				J_P_start_array_list, X12_states);

		// N states in V to D region  transitions
		create_VD_N_part_1_Transitions(vdj, VD_N_states, X1a, X1b, X2, X3,
				V_P_end_states);
		create_VD_N_part_2_Transitions_multi_P(vdj, X3, X5a_states, X5b_states);
		createDTransitions_multi_P(vdj, X5b_states, X6_states, dStates,
				X7a_states, X7b_states, probHolder);

		// N states in D to J region transitions
		create_VD_N_part_1_Transitions_multi_P(vdj, DJ_N_states, X7a_states,
				X7b_states, X8, X9, D_P_end_array_list);
		create_VD_N_part_2_Transitions_multi_P(vdj, X9, X11a_states,
				X11b_states);
		createJTransitions(vdj, X11b_states, X12_states, jStates, probHolder,
				removed_C_Region);

		//*************** set the state transitions *****************//
		// it is important that all the states are added to the model, and that
		// all the state transitions
		// have been created before we attempt to set the transition
		// probabilities
	
		setVTransitions(vdj, vStates, X1a, X1b, probHolder);

		// set the P transitions
		set_P_transitions(vdj, V_P_end_states, probHolder.V_end_P_Fields, X2);

		set_P_transitions_multi_start(vdj, X5a_states, D_P_start_array_list,
				probHolder.D_start_P_Fields, X6_states);
		set_P_transitions_multi_end(vdj, D_P_end_array_list,
				probHolder.D_end_P_Fields, X8);
		set_P_transitions_multi_start(vdj, X11a_states, J_P_start_array_list,
				probHolder.J_start_P_Fields, X12_states);

		// set N transitions - between V and D
		// use probHolder for transition probabilites
		set_VD_N_part_1_Transitions(vdj, VD_N_states, probHolder.VD_N_Fields,
				X1a, X1b, X2, X3, V_P_end_states, probHolder.V_end_P_Fields);

		//set_VD_N_part_2_Transitions(vdj, X3, X5a, X5b_states,
		// D_P_start_states, X6, probHolder.D_start_P_Fields);
		set_VD_N_part_2_Transitions_multi_P(vdj, X3, X5b_states);

		//setDTransitions(vdj, X5b_states, X5a, X6, dStates, X7a, X7b,
		// probHolder);
		setDTransitions_multi_P(vdj, X5b_states, X5a_states, X6_states,
				dStates, X7a_states, X7b_states, probHolder);

		// set N transitions - between D and J
		// names of methods are wrong, because the methods work for both V-D and
		// D-J
		set_VD_N_part_1_Transitions_multi_P(vdj, DJ_N_states,
				probHolder.DJ_N_Fields, X7a_states, X7b_states, X8, X9,
				D_P_end_array_list, probHolder.D_end_P_Fields);

		//set_VD_N_part_2_Transitions(vdj, X9, X11a, X11b_states,
		// J_P_start_states, X12, probHolder.J_start_P_Fields);
		set_VD_N_part_2_Transitions_multi_P(vdj, X9, X11b_states);

		setJTransitions(vdj, X11b_states, X11a_states, X12_states, jStates,
				probHolder, removed_C_Region);

		// close the output Stream
		G_fostream.close();

		// return this completed model to the caller;
		return vdj;

	}//createModel()

	// ******************************************************************* //

	/**
	 * add all the Emission States in "vStates" to the Markov Model "vdj" Note:
	 * states represents a single VGene
	 */
	private void addVStates(SimpleMarkovModel vdj, ArrayList vStates) {
		// add all the states to this model
		try {
			// add the v states
			for (int i = 0; i < vStates.size(); i++) {
				vdj.addState((EmissionState) vStates.get(i));
			}
		} catch (Exception e) {
			throw new Error("Can't add v states to model");
		}//catch
	}// --addVstates()

	// ******************************************************************* //

	/**
	 * add all the Emission States in "dStates" to the Markov Model "vdj" Note:
	 * "dStates" represents multiple genes, so "dStates" is a list of genes,
	 * where each gene is a list of states making up this gene
	 */
	private void addDStates(SimpleMarkovModel vdj, ArrayList dStates) {
		// add all the states in all the genes to this model
		try {

			ArrayList arraylist;

			// for every gene in list
			for (int i = 0; i < dStates.size(); i++) {
				// get current gene
				arraylist = (ArrayList) dStates.get(i);
				// for all the states in current gene
				for (int j = 0; j < arraylist.size(); j++) {
					vdj.addState((EmissionState) arraylist.get(j));
				}
			}
		} catch (Exception e) {
			throw new Error("Can't add D states to model");
		}//catch
	}// --addDstates()

	// ******************************************************************* //

	/**
	 * add all the Emission States in "dStates" to the Markov Model "vdj" Note:
	 * "jStates" represents multiple genes, so "jStates" is a list og genes,
	 * where each gene is a list of states making up this gene
	 */
	private void addJStates(SimpleMarkovModel vdj, ArrayList jStates) {
		// add all the states in all the genes to this model
		try {
			ArrayList arraylist;
			// for every gene in lis
			for (int i = 0; i < jStates.size(); i++) {
				// get current gene
				arraylist = (ArrayList) jStates.get(i);
				// for every state in current gene
				for (int j = 0; j < arraylist.size(); j++) {
					vdj.addState((EmissionState) arraylist.get(j));
				}
			}
		} catch (Exception e) {
			throw new Error("Can't add J states to model");
		}//catch
	}// --addJstates()

	// ******************************************************************* //

	/**
	 * set the emission probability of a N-addition state in the N1-region (V-D
	 * region)
	 */
	private void set_VD_Emission(Distribution nDist) {
		try {
			// set the emission probability for this N state
			nDist.setWeight(DNATools.a(), 0.2360);
			nDist.setWeight(DNATools.c(), 0.2387);
			nDist.setWeight(DNATools.g(), 0.3713);
			nDist.setWeight(DNATools.t(), 0.1538);

		} catch (Exception e) {
			throw new Error("set_VD_Emission: Can't set State Emission");
		}
	}// --set_VD_Emission

	// ******************************************************************* //

	/**
	 * set the emission probability of a N-addition state in the N2-region (D-J
	 * region)
	 */
	private void set_DJ_Emission(Distribution nDist) {
		try {
			// set the emission probability for this N state
			nDist.setWeight(DNATools.a(), 0.2079);
			nDist.setWeight(DNATools.c(), 0.3135);
			nDist.setWeight(DNATools.g(), 0.33);
			nDist.setWeight(DNATools.t(), 0.1485);

		} catch (Exception e) {
			throw new Error("set_DJ_Emission: Can't set State Emission");
		}
	}// --set_DJ_Emission

	
	// *************************** Dot States ************************ //

	/**
	 * Creates an array of DotStates of size "size"
	 */
	private DotState[] createDotStateArray(String name, int size) {
		DotState[] dotStateArray = new DotState[size];

		for (int i = 0; i < size; i++) {

			// create the fake Annotations for the Dot State
			StateInfo stateInfo = new StateInfo("Dot State " + name + i,
					DOT_STATE_TOKEN);
			FakeAnnotation annotation = new FakeAnnotation(stateInfo);

			dotStateArray[i] = new SimpleDotState(name, annotation);
		}

		return dotStateArray;

	}//--createDotStateArray()

	/**
	 * Creates a single DotState with name "name"
	 */
	private DotState createDotState(String name) {
		// create the fake Annotations for the Dot State
		StateInfo stateInfo = new StateInfo("Dot State " + name,
				DOT_STATE_TOKEN);
		FakeAnnotation annotation = new FakeAnnotation(stateInfo);

		return new SimpleDotState(name, annotation);

	}//--createDotStates()

	/**
	 * adds all the Dot States in "dotStateArray" to the Markov Model
	 */
	private void addDotStateArray(MarkovModel vdj, DotState[] dotStateArray) {
		try {
			for (int i = 0; i < dotStateArray.length; i++) {

				vdj.addState(dotStateArray[i]);
			}
		} catch (Exception e) {
			throw new Error("addDotStateArray(): " + e.getMessage());
		}//catch

	}//--addDotStateArray()

	/**
	 * convenience method for adding several pre-defined dot states to the model
	 * at once
	 */
	private void addDotStates(MarkovModel vdj, DotState X1a, DotState X1b,
			DotState X2, DotState X3,
			/*
			 * DotState X5a, DotState X5b, DotState X6, DotState X7a, DotState
			 * X7b,
			 */DotState X8, DotState X9 /*
										 * DotState X11a, DotState X11b,
										 * DotState X12
										 */) {
		try {
			vdj.addState(X1a);
			vdj.addState(X1b);
			vdj.addState(X2);
			vdj.addState(X3);
			vdj.addState(X8);
			vdj.addState(X9);

		} catch (Exception e) {
			throw new Error("addDotStates: " + e.getMessage());
		}//catch
	}// --addDotStates()

	// ******************************************************************* //

	/*
	 * create an aray of p_states belonging to a certain gene (the reverse
	 * compliment of the gene start or end. The emission probabilities are set
	 * so the a P-state can only emit its particular Gene end Reverse Compliment
	 * (RC) nucleotide
	 */
	private ArrayList createPStatesArrayList(int[] advance, FiniteAlphabet dna,
			String stateName, int size, char state_token,
			ArrayList matching_genes_list, boolean gene_start) {
		ArrayList p_results = new ArrayList();
		Sequence gene_seq;
		EmissionState[] p_array;

		// for every gene in gene list, create the RC matching p_array
		for (int i = 0; i < matching_genes_list.size(); i++) {
			gene_seq = (Sequence) matching_genes_list.get(i);

			// create the RC mathcing p array
			p_array = createPStateArray(advance, dna, stateName, size,
					state_token, gene_seq, gene_start);

			// add this p-array to list of p-arrays
			p_results.add(p_array);
		}

		// return the resulting list of p-arrays
		return p_results;

	}//--createPStatesArrayList

	/*
	 * create an aray of p_states belonging to a certain gene (the reverse
	 * compliment of the gene end, either end)
	 */
	private EmissionState[] createPStateArray(int[] advance,
			FiniteAlphabet dna, String stateName, int size, char state_token,
			Sequence matching_gene, boolean gene_start) {
		try {
			Distribution tempDist; // temporary holder
			FakeAnnotation annotation;
			StateInfo stateInfo;
			Symbol sym = null;

			EmissionState[] p_states = new EmissionState[size];

			for (int i = 0; i < size; i++) {

				// create the distribution
				tempDist = DistributionFactory.DEFAULT.createDistribution(dna);
				// add information of this state to the State Info object
				stateInfo = new StateInfo(stateName, state_token);
				annotation = new FakeAnnotation(stateInfo);
				// create the state
				p_states[i] = new SimpleEmissionState(stateName, annotation,
						advance, tempDist);

			}//-- for every p-state we want to create

			//	System.out.println("gene start == " + gene_start);
			// set the emission probability of this array of p-states
			if (gene_start) // if this p_array is at the start of a gene
			{
				//System.out.print("set P_start for: " +
				// matching_gene.getName() + " = ");
				set_P_start(matching_gene, p_states);
			}//--if

			else // this p_array is at the end of a gene
			{
				//System.out.print("set P_end for: " + matching_gene.getName()
				// + " = ");
				set_P_end(matching_gene, p_states);
			}//--else

			// return the array of emission states
			return p_states;

		} catch (Exception e) {
			throw new Error("createPStateArray(): " + e.getMessage());
		}
	}//--createPStateArray

	/**
	 * add the array of P states in "pStates" to the Markov Model "vdj"
	 */
	private void addPStates(MarkovModel vdj, EmissionState[] pStates) {
		try {
			for (int i = 0; i < pStates.length; i++) {
				vdj.addState(pStates[i]);
			}
		} catch (Exception e) {
			throw new Error("addPStates(): " + e.getMessage());
		}//catch
	}// --addPState()

	/**
	 * add the list of arrays of P states in "pStates" to the Markov Model "vdj"
	 */
	private void addPStates_list(MarkovModel vdj, ArrayList pStates_list) {
		try {

			EmissionState[] pState_array;
			for (int i = 0; i < pStates_list.size(); i++) {
				pState_array = (EmissionState[]) pStates_list.get(i);
				for (int j = 0; j < pState_array.length; j++) {
					vdj.addState(pState_array[j]);
				}//--for every pState in array
			}//-- for every array of pStates in ArrayList

		} catch (Exception e) {
			throw new Error("addPStates(): " + e.getMessage());
		}//catch
	}// --addPState()

	/**
	 * create an array of N-addition states of size "size"
	 */
	private EmissionState[] createNStates(int[] advance, FiniteAlphabet dna,
			String geneName, boolean VD, int size) {
		// create all the N states
		try {

			EmissionState[] nStates = new EmissionState[size];

			Distribution tempDist; // temporary holder
			FakeAnnotation annotation;
			StateInfo stateInfo;
			char token;
			for (int i = 0; i < nStates.length; i++) {
				// create the distribution
				tempDist = DistributionFactory.DEFAULT.createDistribution(dna);

				// set the emission probability for this state
				if (VD) // it's a VD-N distribution
				{
					set_VD_Emission(tempDist);
					token = VD_N_TOKEN;
				} else // it's a DJ distribution
				{
					set_DJ_Emission(tempDist);
					token = DJ_N_TOKEN;
				}

				// add information of this state to the State Info object
				stateInfo = new StateInfo(geneName, token);
				annotation = new FakeAnnotation(stateInfo);

				// create the state
				nStates[i] = new SimpleEmissionState(geneName, annotation,
						advance, tempDist);

			}// for

			return nStates;

		} catch (Exception e) {
			throw new Error("createNState(): Can't create N distributions");
		}

	}// --createNStates()

	// ***************************************************************** //

	/**
	 * add the N-addition states in array "n_states" to the MM "vdj"
	 */
	private void addNStates(EmissionState[] n_states, MarkovModel vdj) {
		for (int i = 0; i < n_states.length; i++) {
			try {
				vdj.addState(n_states[i]);
			} catch (Exception ex) {
				throw new Error("addNStates():  " + ex.getMessage());
			}
		}
	}//--addNStates()

	// ***************************************************************** //

	/**
	 * create all the VGene states making up the VGene we are using for our
	 * model (one only)
	 */
	private void createVStates(int[] advance, FiniteAlphabet dna,
			Sequence vSequence, ArrayList vStates,
			ArrayList vStateDistributions, int VGene_start_offset,
			int completeVGeneLength, ProbabilityHolder probHolder,
			double A_probability) {
		// create all the V gene states, and distributions
		try {

			Distribution tempDist; // temporary holder
			FakeAnnotation annotation;
			StateInfo stateInfo;
			Symbol sym = null;

			// when using the VGene nucl position to calculate the probability
			// of mutation,
			// we need to use the absolute nucl. pos in the VGene,
			// in our case, this means adding "i" to the VGene start offset
			int actual_VGene_nucl_pos;

			String seq_string = vSequence.seqString();

			// total probability of there being a mutation from this nucleotide
			double totalMutationProb = probHolder.Gene_Mutation_Field;

			for (int i = 1; i <= vSequence.length(); i++) {

				// create the distribution
				tempDist = DistributionFactory.DEFAULT.createDistribution(dna);
				vStateDistributions.add(tempDist);

				// set the emission probability for this state
				setGeneEmissionProb(vSequence, seq_string, i, tempDist, dna,
						GENE_TYPE_V, completeVGeneLength, VGene_start_offset,
						A_probability);

				// add information of this state to the State Info object
				String geneName = vSequence.getName();
				stateInfo = new StateInfo(geneName, V_STATE_TOKEN, -1,
						(i + VGene_start_offset),
						getChar(vSequence.symbolAt(i)), completeVGeneLength);
				annotation = new FakeAnnotation(stateInfo);

				// create the state
				vStates.add(new SimpleEmissionState(geneName, annotation,
						advance, tempDist));

			}// for
		} catch (Exception e) {
			e.printStackTrace();
			throw new Error("createVStates: Can't create distributions(): "
					+ e.getMessage());
		}
	}// --createVStates()

	// ***************************************************************** //

	/**
	 * create all the DGene states in all the DGenes in "dSequences"
	 */
	private void createDStates(int[] advance, FiniteAlphabet dna,
			ArrayList dSequences, ArrayList dStates,
			ArrayList dStateDistributions, ProbabilityHolder probHolder,
			int completeVGeneLength, double A_probability) {
		Sequence seq; // used as a temporary reference
		String seq_string;
		// create all the D gene states, and distributions

		// probability of there being a mutation at any nucleotide (state) in
		// genes

		double totalMutationProb = probHolder.Gene_Mutation_Field;

		try {

			Distribution tempDist; // temporary holder
			FakeAnnotation annotation;
			StateInfo stateInfo;
			Symbol sym = null;

			// for every DGene
			for (int i = 0; i < dSequences.size(); i++) {
				seq = (Sequence) dSequences.get(i);
				seq_string = seq.seqString(); // string representation of
											  // sequence

				ArrayList stateArrayList = new ArrayList();
				ArrayList distArrayList = new ArrayList();

				// for all the states (nts) in current DGene
				for (int j = 1; j <= seq.length(); j++) {
					// create the distribution
					tempDist = DistributionFactory.DEFAULT
							.createDistribution(dna);
					distArrayList.add(tempDist);

					// set the emission probabilities for this state
					setGeneEmissionProb(seq, seq_string, j, tempDist, dna,
							GENE_TYPE_D, completeVGeneLength, 0, A_probability);

					// add information of this state to the State Info object
					String geneName = seq.getName();
					stateInfo = new StateInfo(geneName, D_STATE_TOKEN, i, j,
							getChar(seq.symbolAt(j)), seq.length());
					annotation = new FakeAnnotation(stateInfo);

					// create the state
					stateArrayList.add(new SimpleEmissionState(geneName,
							annotation, advance, tempDist));
				}
				dStates.add(stateArrayList);
				dStateDistributions.add(distArrayList);
			}// for
		} catch (Exception e) {
			throw new Error("createDStates: Can't create distributions");
		}
	}// --createDStates()

	////////////////////////////////////////////////////////

	/**
	 * get the first character from the name of Symbol "symbol"
	 */
	private char getChar(Symbol symbol) {
		return symbol.getName().charAt(0);
	}//--getChar

	// ***************************************************************** //

	/**
	 * create all the JGene states in all the JGenes in "jSequences"
	 */
	private void createJStates(int[] advance, FiniteAlphabet dna,
			ArrayList jSequences, ArrayList jStates,
			ArrayList jStateDistributions, ProbabilityHolder probHolder,
			int completeVGeneLength, double A_probability) {
		Sequence seq; // used as a temporary reference
		// create all the D gene states, and distributions

		// total probability of there being a mutation from this nucleotide
		double totalMutationProb = probHolder.Gene_Mutation_Field;

		try {

			Distribution tempDist; // temporary holder
			FakeAnnotation annotation;
			StateInfo stateInfo;
			Symbol sym = null;

			// for every JGene
			for (int i = 0; i < jSequences.size(); i++) {
				seq = (Sequence) jSequences.get(i);
				String seq_string = seq.seqString(); // String representation of
													 // the sequence

				ArrayList stateArrayList = new ArrayList();
				ArrayList distArrayList = new ArrayList();

				// for all the states (nts) in current JGene
				for (int j = 1; j <= seq.length(); j++) {
					// create the distribution
					tempDist = DistributionFactory.DEFAULT
							.createDistribution(dna);
					distArrayList.add(tempDist);

					// set the emission probability for this state
					setGeneEmissionProb(seq, seq_string, j, tempDist, dna,
							GENE_TYPE_J, completeVGeneLength, 0, A_probability);

					// add information of this state to the State Info object
					String geneName = seq.getName();
					stateInfo = new StateInfo(geneName, J_STATE_TOKEN, i, j,
							getChar(seq.symbolAt(j)), seq.length());
					annotation = new FakeAnnotation(stateInfo);

					// create the state
					stateArrayList.add(new SimpleEmissionState(geneName,
							annotation, advance, tempDist));
				}
				jStates.add(stateArrayList);
				jStateDistributions.add(distArrayList);
			}// for
		} catch (Exception e) {
			throw new Error("createJStates: Can't create distributions");
		}
	}// --createJStates()

	// ************************************************************** //

	/**
	 * create the State Transitions between the P-States in "pStates" and from
	 * the last P-State for the Dot State "pFinishedState"
	 */
	private void create_P_transitions(SimpleMarkovModel vdj,
			EmissionState[] pStates, DotState pFinishedState) {
		try {
			EmissionState curr_state, next_state;
			next_state = curr_state = null;

			// make transitions from P to P states
			for (int i = 0; i < (pStates.length - 1); i++) // for every p state
														   // but last
			{
				curr_state = pStates[i];
				next_state = pStates[i + 1];

				vdj.createTransition(curr_state, next_state);
				vdj.createTransition(curr_state, pFinishedState);
			}//--for(i)

			// create transition from last p-state to the pFinished state
			vdj.createTransition(next_state, pFinishedState);

		} catch (Exception be) {
			throw new Error("create_P_transitions() : " + be.getMessage());
		}
	}//--create_P_transitions()

	////////////////////////////////////////////////////////////////

	/**
	 * Create the state transition P_start_state -> P1 -> P2 -> PX ->
	 * P_finished_state for every array of P-States in "pStates_list" and their
	 * matching P_start_states and P_finished_states
	 */
	private void create_P_transitions_multi_start(SimpleMarkovModel vdj,
			DotState[] pStartState_list, ArrayList pStates_list,
			DotState[] pFinishedState_list) {
		try {
			EmissionState curr_state, next_state;
			curr_state = next_state = null;
			EmissionState[] pState_array;

			DotState pFinished_state;
			DotState pStart_state;

			// for every array of p-states in array list
			for (int i = 0; i < pStates_list.size(); i++) {
				pState_array = (EmissionState[]) pStates_list.get(i);
				pFinished_state = pFinishedState_list[i]; // get matching end
														  // state for each
														  // p-state
				pStart_state = pStartState_list[i]; // get matching start state
													// for each p-state

				// make transition from X5ai state to X6i state (no exo, no
				// P-add)
				vdj.createTransition(pStart_state, pFinished_state); // no p-add
																	 // and no
																	 // exo if
																	 // this
																	 // path is
																	 // taken

				// make transitions from P to P states in array
				for (int j = 0; j < (pState_array.length - 1); j++) // for every
																	// p state
																	// but first
				{
					curr_state = pState_array[j];
					next_state = pState_array[j + 1];

					vdj.createTransition(pStart_state, curr_state);
					vdj.createTransition(curr_state, next_state);

				}//--for every state in p-state array

				// create transition from start state to last p-state
				vdj.createTransition(pStart_state, next_state);

				// create transition from last p-state to the pFinished state
				vdj.createTransition(next_state, pFinished_state);

			}//--for every p-state-array in array list
		} catch (Exception be) {
			throw new Error("create_P_transitions_multi_start() : "
					+ be.getMessage());
		}
	}//--create_P_transitions_multi_start()

	///////////////////////////////////////////////////////////////

	/**
	 * same as above method, but for P-regions in the End of a Gene
	 */
	private void create_P_transitions_multi_end(SimpleMarkovModel vdj,
			ArrayList pStates_list, DotState X8) {
		try {
			EmissionState curr_state, next_state;
			curr_state = next_state = null;

			EmissionState[] pState_array;

			// for every array of p-states in array list
			for (int i = 0; i < pStates_list.size(); i++) {
				pState_array = (EmissionState[]) pStates_list.get(i);

				// make transitions from P to P states in array
				for (int j = 0; j < (pState_array.length - 1); j++) // for every
																	// p state
																	// but last
				{

					curr_state = pState_array[j];
					next_state = pState_array[j + 1];

					// make transitions from P to P states in array
					vdj.createTransition(curr_state, next_state);

					// make transition from every p_state to X8
					vdj.createTransition(curr_state, X8);

				}//--for every state in p-state array

				// create transition from last p-state to X8
				vdj.createTransition(next_state, X8);

			}//--for every p-state-array in array list

		} catch (Exception be) {
			throw new Error("create_P_transitions_multi_end() : "
					+ be.getMessage());
		}
	}//--create_P_transitions_multi_end()

	// ************************************************************** //

	private void set_P_transitions_multi_start(SimpleMarkovModel vdj,
			DotState[] pStartState_list, ArrayList pStates_list,
			double[] probabilities, DotState[] pFinishedState_list) {
		try {
			EmissionState curr_state, next_state;
			curr_state = next_state = null;
			Distribution dist;

			EmissionState[] pState_array;
			DotState pFinished_state;
			DotState pStart_state;

			double PtoPprob;
			double pFinishedProb;

			// find the summed X5a to Pi prob, so we may normalize the
			// probability
			double X5aToPiSum = getCompleteRelativeProb(probabilities);
			double PaddProb = probabilities[0];

			for (int i = 0; i < pStates_list.size(); i++) {
				pState_array = (EmissionState[]) pStates_list.get(i);
				pFinished_state = pFinishedState_list[i]; // get matching end
														  // state for each
														  // p-state
				pStart_state = pStartState_list[i];

				// set transition from X5ai state to X6i state (no exo, no
				// P-add)
				double pStartTopFinishedProb = 1.0 - probabilities[0];

				dist = vdj.getWeights(pStart_state);
				dist.setWeight(pFinished_state, pStartTopFinishedProb);

				// set transitions from P to P states in array
				for (int j = 0; j < (pState_array.length - 1); j++) // for every
																	// p state
																	// but last
				{

					curr_state = pState_array[j];
					next_state = pState_array[j + 1];

					// set the transition probability of going from X5ai state
					// to Pi state
					double pStartToPiProb = (PaddProb * (getRelativePprob(
							probabilities, j) / X5aToPiSum));
					dist = vdj.getWeights(pStart_state);
					dist.setWeight(curr_state, pStartToPiProb);

					// get the probability of going to next p-addition state
					// form current
					PtoPprob = 1.0; // always

					dist = vdj.getWeights(curr_state);
					dist.setWeight(next_state, PtoPprob);

				}//--for every p-state in array

				// set the transition probability of going from start state to
				// Plast state
				double pStartToPlastProb = (PaddProb * (getRelativePprob(
						probabilities, (probabilities.length - 1)) / X5aToPiSum));

				//		System.out.println(" pStartToPlastProb = " +
				// pStartToPlastProb);

				dist = vdj.getWeights(pStart_state);
				dist.setWeight(next_state, pStartToPlastProb);

				// set transition from last p-state to the pFinished state
				dist = vdj.getWeights(next_state);
				dist.setWeight(pFinished_state, 1.0); // always 1.0

			}//--for every p-state-array in array list
		} catch (Exception be) {
			throw new Error("set_P_transitions_multi_start() : "
					+ be.getMessage());
		}
	}//--set_P_transitions_multi_start()

	/*
	 * helper method to find the probability of going through all states in
	 * "probability array" from offset "offset"
	 */
	private double getRelativePprob(double[] probabilities, int offset) {
		double result = 1.0;
		int reversed_offset = probabilities.length - offset;

		int i;
		for (i = 0; i < reversed_offset; i++) {
			result *= probabilities[i];
		}

		// if there is another possible state, calculate probability of not
		// going to this state
		// and take this into account
		if (i < probabilities.length) {
			double notNextStateProb = 1.0 - probabilities[i];
			result *= notNextStateProb;
		}

		return result;
	}

	/*
	 * the summed probability of all the relative probabilities, so we may
	 * normalize the complete transition
	 */
	private double getCompleteRelativeProb(double[] probabilities) {
		double sum = 0.0;

		for (int i = 0; i < probabilities.length; i++) {
			sum += getRelativePprob(probabilities, i);
		}

		return sum;
	}

	///////////////////////////////////////////////////////////////////////

	private void set_P_transitions_multi_end(SimpleMarkovModel vdj,
			ArrayList pStates_list, double[] probabilities, DotState X8) {
		try {
			EmissionState curr_state, next_state;
			curr_state = next_state = null;
			Distribution dist;

			EmissionState[] pState_array;

			double PtoPprob;
			double pFinishedProb;

			for (int i = 0; i < pStates_list.size(); i++) {
				pState_array = (EmissionState[]) pStates_list.get(i);

				// set transitions from P to P states in array
				for (int j = 0; j < (pState_array.length - 1); j++) // for every
																	// p state
																	// but last
				{

					// set the probabilities of transitions from P to P states

					curr_state = pState_array[j];
					next_state = pState_array[j + 1];

					dist = vdj.getWeights(curr_state);

					PtoPprob = probabilities[j + 1];
					pFinishedProb = 1.0 - PtoPprob; // probability of end of p
													// addition

					dist.setWeight(next_state, PtoPprob);
					dist.setWeight(X8, pFinishedProb);

				}//--for every p-state in array

				// set transition from last p-state to the pFinished state
				dist = vdj.getWeights(next_state);
				dist.setWeight(X8, 1.0); // always 1.0

			}//--for every p-state-array in array list
		} catch (Exception be) {
			throw new Error("set_P_transitions_multi_end() : "
					+ be.getMessage());
		}
	}//--set_P_transitions_multi_end()

	// ************************************************************** //

	private void set_P_transitions(SimpleMarkovModel vdj,
			EmissionState[] pStates, double[] probabilities,
			DotState pFinishedState) {
		try {
			EmissionState curr_state, next_state;
			Distribution dist;

			double PtoPprob;
			double pFinishedProb;

			// set the probabilities of transitions from P to P states
			for (int i = 0; i < (pStates.length - 1); i++) // for every p state
														   // but last
			{
				curr_state = pStates[i];
				next_state = pStates[i + 1];

				dist = vdj.getWeights(curr_state);

				PtoPprob = probabilities[i + 1];
				pFinishedProb = 1.0 - PtoPprob; // probability of end of p
												// addition

				dist.setWeight(next_state, PtoPprob);
				dist.setWeight(pFinishedState, pFinishedProb);

			}//--for(i)
		} catch (Exception be) {
			throw new Error("set_P_transitions() : " + be.getMessage());
		}
	}//--set_P_transitions()

	// ************************************************************** //
	// ********************** create N transitions ******************** //
	// ************************************************************** //

	/**
	 * create the transitions surrounding the VD-N region (first part), assuming
	 * there are multiple P-regions (one for each gene the region connects to)
	 */
	private void create_VD_N_part_1_Transitions_multi_P(SimpleMarkovModel vdj,
			EmissionState[] VD_N_53_states, DotState[] X7a_states,
			DotState[] X7b_states, DotState X8, DotState X9,
			ArrayList P_D_end_states_array_list) {
		try {

			// local vars
			EmissionState currN = null;
			EmissionState nextN = null;

			EmissionState[] P_D_end_state_array;

			// make the X7ai to first P_Vi state and X7ai to X8 transitions

			// for every P_D_end state array in array list
			for (int i = 0; i < P_D_end_states_array_list.size(); i++) {
				// get state array
				P_D_end_state_array = (EmissionState[]) P_D_end_states_array_list
						.get(i);

				// make the X7ai to first state in P_V_end state array
				vdj.createTransition(X7a_states[i], P_D_end_state_array[0]);

				// make the X7ai to X8 transitions
				vdj.createTransition(X7a_states[i], X8);

				// make the X7bi to X8 transition
				vdj.createTransition(X7b_states[i], X8);
			}

			// already made Pi to X8 and P-last to X8 transitions in create P
			// transitions

			// make the X8 to N1 and X8 to X9 transition
			vdj.createTransition(X8, X9);
			vdj.createTransition(X8, VD_N_53_states[0]);

			// make the Nn to Nn+1 transitions, and Nn to X9 transitions
			for (int i = 0; i < (VD_N_53_states.length - 1); i++) // for all but
																  // the last N
			{
				currN = VD_N_53_states[i];
				nextN = VD_N_53_states[i + 1];

				vdj.createTransition(currN, nextN);

				// make transitions to X9
				vdj.createTransition(currN, X9);
			}

			// make transitions from last N to X9
			vdj.createTransition(nextN, X9); // nextN = last N

		} catch (Exception e) {
			throw new Error("create_VD_N_part_1_Transitions_multi_P():  "
					+ e.getMessage());
		}//catch
	}// --create_VD_N_part_1_Transitions_multi_P

	///////////////////////////////////////////////////////////////////

	/**
	 * create the state transitions surronding the VD-N region (assuming only
	 * one P-region for each of the four P-region (not one for each gene))
	 */
	private void create_VD_N_part_1_Transitions(SimpleMarkovModel vdj,
			EmissionState[] VD_N_53_states, DotState X1a, DotState X1b,
			DotState X2, DotState X3, EmissionState[] P_V_end_states) {
		try {

			// local vars
			EmissionState currN = null;
			EmissionState nextN = null;

			// make the X1a to first P_V state and X1a to X2 transitions
			vdj.createTransition(X1a, P_V_end_states[0]);
			vdj.createTransition(X1a, X2);

			// already made P-last to X2 in create P transitions

			// make the last P_V state to X2 transition
			//int last_p_index = P_V_end_states.length - 1;
			//vdj.createTransition(P_V_end_states[last_p_index], X2);

			// make the X1b to X2 transition
			vdj.createTransition(X1b, X2);

			// make the X2 to N1 and X2 to X3 transition
			vdj.createTransition(X2, X3);
			vdj.createTransition(X2, VD_N_53_states[0]);

			// make the Nn to Nn+1 transitions, and Nn to X3 transitions
			for (int i = 0; i < (VD_N_53_states.length - 1); i++) // for all but
																  // the last N
			{
				currN = VD_N_53_states[i];
				nextN = VD_N_53_states[i + 1];

				vdj.createTransition(currN, nextN);

				// make transitions to X3
				vdj.createTransition(currN, X3);
			}

			// make transitions from last N to X3
			vdj.createTransition(nextN, X3); // nextN = last N

		} catch (Exception e) {
			throw new Error("create_VD_N_part_1_Transitions():  "
					+ e.getMessage());
		}//catch
	}// --create_VD_N_part_1_Transitions

	// ************************************************************** //

	/**
	 * create the transitions surrounding the VD region (part 2) when we have
	 * multiple P-addition regions (one for each corresponding gene)
	 */
	private void create_VD_N_part_2_Transitions_multi_P(SimpleMarkovModel vdj,
			DotState X3, DotState[] X5a_states, DotState[] X5b_states) {
		try {

			// make the X3 to X5b's and
			// make transition from X5b states to X5a states
			for (int i = 0; i < X5b_states.length; i++) {
				vdj.createTransition(X3, X5b_states[i]); // always taken
				vdj.createTransition(X5b_states[i], X5a_states[i]); // no exo
																	// and maybe
																	// p-add if
																	// this path
																	// is taken
			}

			// already made X5a's to X6's, X5aToPi's, PitoPi+1 and P-last to X6
			// in create_P_transitions()

		} catch (Exception e) {
			throw new Error(
					"create_VD_N_part_2_Transitions_multi_P(): Can't create transitions");
		}//catch
	}// --create_VD_N_part_2_Transitions_multi_P()

	////////////////////////////////////////////////////////////////////////

	/**
	 * create the transitions surrounding the VD-N region when there is only one
	 * P-addition region (not multiple P regions for each corresponding gene)
	 */
	private void create_VD_N_part_2_Transitions(SimpleMarkovModel vdj,
			DotState X3, DotState X5a, DotState[] X5b_states, DotState X6,
			EmissionState[] P_D_states) {
		try {

			// make the X3 to X5b's and
			// make transition from every X5b state to the X5a state
			for (int i = 0; i < X5b_states.length; i++) {
				vdj.createTransition(X3, X5b_states[i]); // always taken
				vdj.createTransition(X5b_states[i], X5a); // no exo and maybe
														  // p-add if this path
														  // is taken
			}

			// make the X5a to first P_D state
			vdj.createTransition(X5a, P_D_states[0]);

			vdj.createTransition(X5a, X6); // no p-add and no exo if this path
										   // is taken

			// already made PtoP and P-last to X6 in create_P_transitions()

		} catch (Exception e) {
			throw new Error("Can't create transitions");
		}//catch
	}// --create_VD_N_part_2_Transitions()

	// ************************************************************** //
	// ********************** Set N transitions ******************** //
	// ************************************************************** //

	/**
	 * set the probabilities of the state transitions created by
	 * "create_VD_N_part_1_Transitions_multi_P(..)"
	 */
	private void set_VD_N_part_1_Transitions_multi_P(SimpleMarkovModel vdj,
			EmissionState[] N_53_states, double[] N_53_Probs,
			DotState[] X7a_states, DotState[] X7b_states, DotState X8,
			DotState X9, ArrayList P_end_states_list, double[] P_probs) {
		try {

			//local vars
			EmissionState currN = null;
			EmissionState nextN = null;
			Distribution dist;

			EmissionState[] P_end_state_array;

			// set the X7ai to start of P_end's and X7ai to X8 transitions

			// these probabilities are valid for the entire for loop below

			double X7aiToP_endProb = P_probs[0]; // probability of p-add given
												 // there is no exo. act.
			double X7aiToX8Prob = 1.0 - X7aiToP_endProb; // dependant

			// set the X7bi to X8 transition
			double X7biToX8Prob = 1.0; // always 1.0

			// for every P_end state array in arrayList
			for (int i = 0; i < P_end_states_list.size(); i++) {
				// get p_state array
				P_end_state_array = (EmissionState[]) P_end_states_list.get(i);

				// set X7ai to Pi[0] transition and X7ai to X8 transition
				dist = vdj.getWeights(X7a_states[i]);
				dist.setWeight(P_end_state_array[0], X7aiToP_endProb);
				dist.setWeight(X8, X7aiToX8Prob);

				// set X7bi to X8 transition
				dist = vdj.getWeights(X7b_states[i]);
				dist.setWeight(X8, X7biToX8Prob);
			}

			// set the X8 to N1 and X8 to X9 transition
			double X8ToN1Prob = N_53_Probs[0];
			double X8ToX9Prob = 1.0 - X8ToN1Prob; // dependant, X8 to X9 skips
												  // the entire N addition

			dist = vdj.getWeights(X8);
			dist.setWeight(N_53_states[0], X8ToN1Prob);
			dist.setWeight(X9, X8ToX9Prob);

			// set the Nn to Nn+1 transitions, and Nn to X9 transitions

			double NToNProb;
			double NToX9Prob;

			for (int i = 0; i < (N_53_states.length - 1); i++) // for every
															   // state but the
															   // last
			{
				currN = N_53_states[i];
				nextN = N_53_states[i + 1];

				NToNProb = N_53_Probs[i + 1];
				NToX9Prob = 1.0 - NToNProb; // dependant

				dist = vdj.getWeights(currN);
				dist.setWeight(nextN, NToNProb);
				dist.setWeight(X9, NToX9Prob);
			}

			// set the lastN to X9 transition
			dist = vdj.getWeights(nextN); // last N
			dist.setWeight(X9, 1.0); // always 1.0

		} catch (Exception e) {
			throw new Error("set_VD_N_part_1_Transitions_multi_P():  "
					+ e.getMessage());
		}//catch
	}//--set_VD_N_part_1_Transitions_multi_P

	////////////////////////////////////////////////////////////////

	/**
	 * set the probabilities of the state transitions created by
	 * "create_VD_N_part_1_Transitions(...)"
	 */
	private void set_VD_N_part_1_Transitions(SimpleMarkovModel vdj,
			EmissionState[] N_53_states, double[] N_53_Probs, DotState X1a,
			DotState X1b, DotState X2, DotState X3,
			EmissionState[] P_end_states, double[] P_probs) {
		try {

			//local vars
			EmissionState currN = null;
			EmissionState nextN = null;
			Distribution dist;

			// set the X1a to start of P_end and X1a to X2 transition
			//double X1aToP_endProb = 0.5;
			double X1aToP_endProb = P_probs[0]; // probability of p-add given
												// there is no exo. act.

			double X1aToX2Prob = 1.0 - X1aToP_endProb; // dependant

			dist = vdj.getWeights(X1a);
			dist.setWeight(P_end_states[0], X1aToP_endProb);
			dist.setWeight(X2, X1aToX2Prob);

			// already set in set_P_Transitions()

			// set the P_V_end to X2 transition
			//double P_endToX2Prob = 1.0 - PtoPProb;
			//dist = vdj.getWeights(P_end);
			//dist.setWeight(P_end, PtoPProb);
			//dist.setWeight(X2, P_endToX2Prob);

			// set the X1b to X2 transition
			double X1bToX2Prob = 1.0; // always 1.0

			dist = vdj.getWeights(X1b);
			dist.setWeight(X2, X1bToX2Prob);

			// set the X2 to N1 and X2 to X3 transition
			double X2ToN1Prob = N_53_Probs[0];
			double X2ToX3Prob = 1.0 - X2ToN1Prob; // dependant, X2 to X3 skips
												  // the entire N addition

			dist = vdj.getWeights(X2);
			dist.setWeight(N_53_states[0], X2ToN1Prob);
			dist.setWeight(X3, X2ToX3Prob);

			// set the Nn to Nn+1 transitions, and Nn to X3 transitions

			double NToNProb;
			double NToX3Prob;

			for (int i = 0; i < (N_53_states.length - 1); i++) // for every
															   // state but the
															   // last
			{
				currN = N_53_states[i];
				nextN = N_53_states[i + 1];

				NToNProb = N_53_Probs[i + 1];
				NToX3Prob = 1.0 - NToNProb; // dependant

				dist = vdj.getWeights(currN);
				dist.setWeight(nextN, NToNProb);
				dist.setWeight(X3, NToX3Prob);
			}

			// set the lastN to X3 transition
			dist = vdj.getWeights(nextN); // last N
			dist.setWeight(X3, 1.0); // always 1.0

		} catch (Exception e) {
			throw new Error("set_VD_N_part_1_Transitions():  " + e.getMessage());
		}//catch
	}//--set_VD_N_part_1_Transitions

	// ************************************************************** //

	/**
	 * set the probabilities of the transitions created by
	 * "create_VD_N_part_2_Transitions_multi_P(...)'
	 */
	private void set_VD_N_part_2_Transitions_multi_P(SimpleMarkovModel vdj,
			DotState X3, DotState[] X5b_states) {
		try {

			//local vars
			EmissionState currN = null;
			EmissionState nextN = null;
			Distribution dist;

			// the length of the "X5b_states" array is equal to the number of
			// dGenes (1 to 1 relationship)
			final double NUMBER_OF_DGENES = X5b_states.length;

			// when setting the transition probabilities from a DotState into a
			// start
			// state of a DGene, we must take into account whether or not the
			// Gene
			// is inversed, as inverted DGenes are less likely than non-inverted
			// DGenes
			// also, we know that every DGene has an inversion, so half the
			// "dGeneCount" is
			// inverted DGenes. Also, every second DGene is inverted, so we
			// don't have to test
			// each DGene to see if it's inverted
			final double INVERTED_DGENE_PROB = (double) 9 / (double) 246; // the
																		  // probability
																		  // of
																		  // an
																		  // inverted
																		  // DGene
																		  // being
																		  // used
																		  // in a
																		  // sequence
			final double NON_INVERTED_DGENE_PROB = (double) 1.0
					- INVERTED_DGENE_PROB; // probability of a non-inverted

			final double NORMALIZED_NON_INVERTED_PROB = NON_INVERTED_DGENE_PROB
					/ NUMBER_OF_DGENES; // the normalized non-inverted
										// probability for a transition
			final double NORMALIZED_INVERTED_PROB = INVERTED_DGENE_PROB
					/ NUMBER_OF_DGENES; // the normalized inverted probability
										// for a transition

			/*
			 * System.out.println("number of DGenes: " + NUMBER_OF_DGENES);
			 * System.out.println("inverted normalized prob: " +
			 * NORMALIZED_INVERTED_PROB); System.out.println("non inverted
			 * normalized prob: " + NORMALIZED_NON_INVERTED_PROB);
			 * 
			 * BufferedReader streami = new BufferedReader( new
			 * InputStreamReader(System.in));
			 * 
			 * streami.readLine();
			 */

			// set the X3 to every X5b state transition
			// the length of the "X5b_states" array is equal to the number of
			// dGenes (1 to 1 relationship)
			double X3ToX5bProb = 1.0 / X5b_states.length; // constant (1.0 /
														  // number of states)
			dist = vdj.getWeights(X3);

			boolean invertedDGene = false; // first DGene is not inverted, but
										   // every second one is

			// for every X5b state
			for (int i = 0; i < X5b_states.length; i++) {
				// dist.setWeight(X5b_states[i], X3ToX5bProb);

				// test if current transition (effectively into a DGene) should
				// have the inverted or non-inverted probability
				if (invertedDGene) {
					dist.setWeight(X5b_states[i], NORMALIZED_INVERTED_PROB);
					invertedDGene = false; // change the inverted dgene flag
				} else {
					dist.setWeight(X5b_states[i], NORMALIZED_NON_INVERTED_PROB);
					invertedDGene = true; // change the inverted dgene flag
				}
			}//--for(i)

			//*** X5b's to X5a's transition probability is set in
			// setDtransitions() and setVtransitions()

			// X5a's to X6's, X5a's to Pi's, PtoP and P-last to X6 is already
			// set in set_P_transitions()

		} catch (Exception e) {
			throw new Error("set_VD_N_part_2_Transitions_multi_P():  "
					+ e.getMessage());
		}//catch
	}// --set_VD_N_part_2_Transitions_multi_P()

	/////////////////////////////////////////////////////////////////

	/**
	 * set the probabilities for the transitions created by
	 * "create_VD_N_part_2_Transitions()",
	 */
	private void set_VD_N_part_2_Transitions(SimpleMarkovModel vdj,
			DotState X3, DotState X5a, DotState[] X5b_states,
			EmissionState[] P_start_states, DotState X6, double[] P_probs) {
		try {

			//local vars
			EmissionState currN = null;
			EmissionState nextN = null;
			Distribution dist;

			// set the X3 to every X5b state transition
			double X3ToX5bProb = 1.0 / X5b_states.length; // constant (1.0 /
														  // number of states)
			dist = vdj.getWeights(X3);

			// for every X5b state
			for (int i = 0; i < X5b_states.length; i++) {
				dist.setWeight(X5b_states[i], X3ToX5bProb);
			}

			//*** X5b to X5a transition probability is set in setDtransitions()
			// and setVtransitions()

			// set the X5a to P_start, and X5a to X6
			double X5aToP_startProb = P_probs[0]; // probability of p-addition
												  // given there is no exo. act.
			double X5aToX6Prob = 1.0 - X5aToP_startProb; // prob of no p-add,
														 // given no exo act.

			dist = vdj.getWeights(X5a);
			dist.setWeight(X6, X5aToX6Prob);
			dist.setWeight(P_start_states[0], X5aToP_startProb);

			// PtoP and P-last to X6 is already set in set_P_transitions()

		} catch (Exception e) {
			throw new Error("set_VD_N_part_2_Transitions():  " + e.getMessage());
		}//catch
	}// --set_VD_N_part_2_Transitions

	// ************************************************************** //
	// ********************* V transitions ************************* //
	// ************************************************************** //

	/**
	 * create the transitions between the V states and their surrounding dot
	 * states
	 */
	private void createVTransitions(SimpleMarkovModel vdj, ArrayList vStates,
			DotState X1a, DotState X1b, ProbabilityHolder probHolder) {
		try {
			// make Magical to V0, V to V, Vlast to X1a and Vn to X1b
			EmissionState currVState = null; // must be initialized
			EmissionState nextVState = null;

			// Magical to V0
			vdj.createTransition(vdj.magicalState(), (EmissionState) vStates
					.get(0));

			// first find the family of gene
			EmissionState vState = (EmissionState) vStates.get(0);
			Annotation anno = vState.getAnnotation();
			StateInfo si = (StateInfo) anno.getProperty(null);
			int familyIndex = findFamilyIndex(si.geneName);

			// then create the probability array from the mean and std deviation
			double mean = probHolder.V_end_exo_mean_Fields[familyIndex - 1];
			double stdDev = probHolder.V_end_exo_stdDev_Fields[familyIndex - 1];
			double[] V_end_exo_probs = probHolder.getExoProbArray(mean, stdDev,
					MIN_PROB_LIMIT, true);

			// create transitions V to V and Vx to X1b where there is exo nucl.
			// activity
			Distribution dist;
			for (int i = 0; i < (vStates.size() - 1); i++) // for every V state
														   // but last
			{
				currVState = (EmissionState) vStates.get(i);

				// V to V transition
				nextVState = (EmissionState) vStates.get(i + 1);
				vdj.createTransition(currVState, nextVState);
				// Model exo nucl. actitvity
				if ((vStates.size() - i) <= V_end_exo_probs.length) {
					vdj.createTransition(currVState, X1b);
				}
			}//--for(i)
			// make connections from last V state
			vdj.createTransition(nextVState, X1a);

		} catch (Exception e) {
			throw new Error("createVTransitions(): " + e.getMessage());
		}//catch
	}// --createVTransitions()

	////////////////////////////////////////////////////////////////////

	/**
	 * get the gene family of a Gene with name "geneName" by parsing the name
	 * returns an integer between 1 and 7 representing the gene Name (assuming
	 * nothing went wrong)
	 */
	private int findFamilyIndex(String geneName) {
		// all d-genes except for the DIR D-genes start with IGH
		if (geneName.indexOf("DIR") != -1) // is this a DIR D-gene
		{
			// use the IGHD 5 family for the DIR D-Genes as they do not have
			// their own unique probabilities
			return 5;
		} else // should be IGHD gene
		{
			int startIGH = geneName.indexOf("IGH"); // find the location of the
													// gene family number
			if (startIGH == -1) // no "IGH" string present
			{
				throw new Error(
						"Invalid D-Gene sequence name:  neiher IGH nor DIR");
			}
			int familyNumberIndex = startIGH + 4;
			char number = geneName.charAt(familyNumberIndex);

			// parse the character where the gene family number should be
			// located
			int result = Integer.parseInt("" + number);

			// return gene family number
			return result;
		}

	}//--findFamilyIndex()

	////////////////////////////////////////////////////////////////////

	/**
	 * set the probabilities of the transitions created in createVTransitions
	 */
	private void setVTransitions(SimpleMarkovModel vdj, ArrayList vStates,
			DotState X1a, DotState X1b, ProbabilityHolder probHolder) {
		try {
			// find the exo. nucl. probability based on the family of gene
			EmissionState vState = (EmissionState) vStates.get(0);
			Annotation anno = vState.getAnnotation();
			StateInfo si = (StateInfo) anno.getProperty(null);
			int familyIndex = findFamilyIndex(si.geneName);

			//System.out.println("family Index = " + familyIndex);

			double mean = probHolder.V_end_exo_mean_Fields[familyIndex - 1];

			//System.out.println("mean = " + mean);

			double stdDev = probHolder.V_end_exo_stdDev_Fields[familyIndex - 1];

			//double [] V_end_exo_probs = probHolder.getExoProbArray(mean,
			// stdDev, MIN_PROB_LIMIT, true);
			double[] V_end_exo_probs = probHolder.getExoProbArray(mean, stdDev,
					MIN_PROB_LIMIT, false);

			// *** Model exonuclease activity by Setting connections/transitions
			// from the last X States (except last) in the V Gene to X1b
			EmissionState currVState = null;
			EmissionState nextVState = null;
			Distribution dist;

			// Set the Magical to V0 transition Probability
			double MagicalToVProb = 1.0;
			dist = vdj.getWeights(vdj.magicalState());
			dist.setWeight((EmissionState) vStates.get(0), MagicalToVProb);

			double temp;
			double VToVProb; // 1.0 - prob of exo nuclease activity from state

			int v_size = vStates.size(); // calculate it once only

			// used to traverse the exo probability array
			//( first probability is for 0 removals ), so start at 1
			int probArrayIndex = 1;

			for (int i = 0; i < (v_size - 1); i++) // for every V state but last
			{
				currVState = (EmissionState) vStates.get(i);
				nextVState = (EmissionState) vStates.get(i + 1);
				dist = vdj.getWeights(currVState);

				/*
				 * // test for exo. nucl. at end if((v_size - i) <=
				 * V_end_exo_probs.length) // there is a exo nucl. transition {
				 * temp = V_end_exo_probs[probArrayIndex++]; VToVProb = 1.0 -
				 * temp;
				 * 
				 * dist.setWeight(nextVState, VToVProb);
				 * dist.setWeight(X1b,temp); System.out.println("V to X1b Prob: " +
				 * i + " = " + temp); // System.out.println("V to V Prob" + i + "
				 * to = " + VToVProb); } else { // set the V to V probability
				 * where there is no exo transition dist.setWeight(nextVState,
				 * 1.0); }
				 */

				// test for exo. nucl. at end
				if ((v_size - i) <= V_end_exo_probs.length) // there is a exo
															// nucl. transition
				{
					temp = V_end_exo_probs[v_size - i - 1];
					VToVProb = 1.0 - temp;

					dist.setWeight(nextVState, VToVProb);
					dist.setWeight(X1b, temp);
					//	 System.out.println("V to X1b Prob: " + i + " = " + temp);
					// System.out.println("V to V Prob" + i + " to = " +
					// VToVProb);
				} else {
					// set the V to V probability where there is no exo
					// transition
					dist.setWeight(nextVState, 1.0);
				}
			}// --for(i)

			// make connections from last V state to X1a
			dist = vdj.getWeights(nextVState);
			dist.setWeight(X1a, 1.0);

		} catch (Exception e) {
			throw new Error("SetVTransitions(): " + e.getMessage());
		}//catch
	}// --setVTransitions

	// ************************************************************** //
	// ***************** D transitions ********************** //
	// ************************************************************** //

	/**
	 * create the transitions between the D states and their surrounding dot
	 * states when we have separate P-state-regions for each D Gene
	 */
	private void createDTransitions_multi_P(SimpleMarkovModel vdj,
			DotState[] X5b_states, DotState[] X6_states, ArrayList dStates,
			DotState[] X7a_states, DotState[] X7b_states,
			ProbabilityHolder probHolder) {
		try {
			// local vars
			ArrayList dSequence; // one gene sequence
			EmissionState dState;

			// array holding the exo nucl. probabilities
			double[] D_start_exo_probs = null;
			double[] D_end_exo_probs = null;

			Annotation anno;
			StateInfo si;
			int familyIndex;
			double mean_start, stdDev_start, mean_end, stdDev_end;

			// make the X6i to DX1's and X5bi to DXn+1 transitions
			// make the DXlast to X7ai and DXnx to X7bi transition
			int dSize = 0;
			for (int i = 0; i < dStates.size(); i++) {
				// get the current gene and the first state of this gene
				dSequence = (ArrayList) dStates.get(i);
				dSize = dSequence.size();
				// get first state of d-gene
				dState = (EmissionState) dSequence.get(0); // get DX1 state

				// find the exo. nucl. prob. arrays based on family of Gene
				anno = dState.getAnnotation();
				si = (StateInfo) anno.getProperty(null);
				familyIndex = findFamilyIndex(si.geneName);

				// exo start array
				mean_start = probHolder.D_start_exo_mean_Fields[familyIndex - 1];
				stdDev_start = probHolder.D_start_exo_stdDev_Fields[familyIndex - 1];
				D_start_exo_probs = probHolder.getExoProbArray(mean_start,
						stdDev_start, MIN_PROB_LIMIT, false);

				// exo end array
				mean_end = probHolder.D_end_exo_mean_Fields[familyIndex - 1];
				stdDev_end = probHolder.D_end_exo_stdDev_Fields[familyIndex - 1];
				D_end_exo_probs = probHolder.getExoProbArray(mean_end,
						stdDev_end, MIN_PROB_LIMIT, true);

				// X6i to DX1 transitions and X5bi to DXn+1 transitions
				vdj.createTransition(X6_states[i], dState);

				for (int j = 1; j < (dSize - 1); j++) // for every d state but
													  // first and last
				{
					dState = (EmissionState) dSequence.get(j);

					// test for start exo transition
					if (j < D_start_exo_probs.length)
						vdj.createTransition(X5b_states[i], dState);

					// test for end exo transition
					if (dSize - j <= D_end_exo_probs.length)
						vdj.createTransition(dState, X7b_states[i]);
				}//--for(j)
			}//--for(i)

			// *** create D to D transitions *** //

			EmissionState currD;
			EmissionState nextD = null;

			for (int i = 0; i < dStates.size(); i++) {
				// create transitions
				dSequence = (ArrayList) dStates.get(i);

				for (int j = 0; j < (dSequence.size() - 1); j++) // for every D
																 // state but
																 // last
				{
					currD = (EmissionState) dSequence.get(j);
					nextD = (EmissionState) dSequence.get(j + 1);

					vdj.createTransition(currD, nextD);
				}//--for(j)

				// create DXlast's to X7ai transitions
				vdj.createTransition(nextD, X7a_states[i]); // nextD equals last
															// D state
			}//--for(i)

		} catch (Exception e) {
			throw new Error("createDTransitions_multi_P(): " + e.getMessage());
		}//catch
	}// --createDTransitions_multi_P()



	// ************************************************************** //

	/**
	 * set the probabilities of the transitions created in
	 * createDTransitions_multi_P(..)
	 */
	private void setDTransitions_multi_P(SimpleMarkovModel vdj,
			DotState[] X5b_states, DotState[] X5a_states, DotState[] X6_states,
			ArrayList dStates, DotState[] X7a_states, DotState[] X7b_states,
			ProbabilityHolder probHolder) {
		try {

			// find the exo. nucl. probability based on the family of gene
			double[] D_start_exo_probs = null;
			double[] D_end_exo_probs = null;

			Annotation anno;
			StateInfo si;
			int familyIndex;

			double mean_start, mean_end, stdDev_start, stdDev_end;

			// local vars
			ArrayList dSequence; // one gene sequence
			Distribution dist;
			EmissionState dState;
			EmissionState currD = null;
			EmissionState nextD = null;

			// set the X6's to DX1's, X5b's to X5as and X5b's DXn+1 transition,
			// plus Dx to Dx+1 probs

			double exo_prob = 0; // used to calculate the remainding transition
								 // prob from currD to nextD
			// when there is exo activity involved
			double remainder = 0; // the remainder used for DX to DX+1 trans.
								  // probs.

			int dSize; // size of dSequence
			int D_end_exo_prob_index = 1; // used to keep track of exo nucl.
										  // act. probs at end of D genes

			int dGeneCount = dStates.size(); // total number of dGenes (
											 // important for calculating
											 // probabilities

			for (int i = 0; i < dStates.size(); i++) {
				dSequence = (ArrayList) dStates.get(i);
				dSize = dSequence.size();
				dState = (EmissionState) dSequence.get(0);

				// find the exo. nucl. prob. based on family of Gene
				anno = dState.getAnnotation();
				si = (StateInfo) anno.getProperty(null);
				familyIndex = findFamilyIndex(si.geneName);

				// exo start
				mean_start = probHolder.D_start_exo_mean_Fields[familyIndex - 1];
				stdDev_start = probHolder.D_start_exo_stdDev_Fields[familyIndex - 1];
				D_start_exo_probs = probHolder.getExoProbArray(mean_start,
						stdDev_start, MIN_PROB_LIMIT, false);

				// exo end
				mean_end = probHolder.D_end_exo_mean_Fields[familyIndex - 1];
				stdDev_end = probHolder.D_end_exo_stdDev_Fields[familyIndex - 1];
				//D_end_exo_probs = probHolder.getExoProbArray(mean_end,
				// stdDev_end, MIN_PROB_LIMIT, true);
				D_end_exo_probs = probHolder.getExoProbArray(mean_end,
						stdDev_end, MIN_PROB_LIMIT, false);

				// find probability of no exo nucl. activity at the start and at
				// the end of gene
				double no_exo_start_prob = D_start_exo_probs[0]; // the
																 // probability
																 // of zero
																 // removals

				//System.out.println("no exo prob = " + no_exo_start_prob);

				double X5biToX5aiProb = no_exo_start_prob; // 1 to 1
														   // relationship
				double X6itoDX1sProb = 1.0; // 1 to 1 relationship prob

				// set the X5b's to X5a's (the probability of noe exo nucl.
				// actitivity)
				dist = vdj.getWeights(X5b_states[i]); // curr X5b state
				dist.setWeight(X5a_states[i], X5biToX5aiProb);

				// set X6i to DX1i transition prob
				dist = vdj.getWeights(X6_states[i]);
				dist.setWeight(dState, X6itoDX1sProb); // use the divided
													   // probability

				// set DX1 to DX2
				dist = vdj.getWeights(dState);
				dist.setWeight((EmissionState) dSequence.get(1), 1.0);

				D_end_exo_prob_index = 1; // reset for every gene

				for (int j = 1; j < (dSequence.size() - 1); j++) // for every D
																 // state but
																 // first and
																 // last
				{
					currD = (EmissionState) dSequence.get(j);
					nextD = (EmissionState) dSequence.get(j + 1);

					remainder = 1.0; // the complete probability

					// test for exo nuclease activity at start and end, they may
					// be overlapping

					// exo start
					if (j < D_start_exo_probs.length) { // we have exo
														// actitivity at the
														// start
						exo_prob = D_start_exo_probs[j]; // start at 1, as 1 = 1
														 // removal
						dist = vdj.getWeights(X5b_states[i]); // curr X5b state
						dist.setWeight(currD, exo_prob); // 1 to 1 relationship
					}

					// exo end
					if ((dSize - j) <= D_end_exo_probs.length) { // we have exo
																 // actitivity
																 // at the end
						exo_prob = D_end_exo_probs[(dSize - j - 1)]; // index
																	 // into
																	 // reversed
																	 // array
						remainder -= exo_prob; // the remainder of the
											   // probability is prob. of going
											   // to the next D state

						dist = vdj.getWeights(currD);
						dist.setWeight(X7b_states[i], exo_prob);
					}

					// set the D to next D transition
					dist = vdj.getWeights(currD);
					dist.setWeight(nextD, remainder);

				}//--for state in gene

				// set the DXlast to X7a transition
				dist = vdj.getWeights(nextD);
				dist.setWeight(X7a_states[i], 1.0); // always 1.0

			}//--for genes in arraylist

		} catch (Exception e) {
			throw new Error("setDTransitions_multi_P():  " + e.getMessage());
		}//catch
	}// --setDTransitions_multi_P()

	

	// ************************************************************** //
	// ************** J Transitions ******************* //
	// ************************************************************** //

	/**
	 * create the transitions between the J states and their surrounding dot
	 * states, some of the transitions depends on the variable
	 * "removed_C_Region", which allows for missing parts of the JGene
	 */
	private void createJTransitions(SimpleMarkovModel vdj,
			DotState[] X11b_states, DotState[] X12_states, ArrayList jStates,
			ProbabilityHolder probHolder, boolean removed_C_Region) //,
																	// EmissionState
																	// JC_bridge_state)
	{
		try {
			// local vars
			ArrayList jSequence; // one gene sequence
			EmissionState jState;
			EmissionState currJ = null;
			EmissionState nextJ = null;

			double[] J_start_exo_probs = null;
			Annotation anno;
			StateInfo si;
			int familyIndex;
			double mean_start, stdDev_start;

			// make the X12's to JX1's and X11b's to JXn+1 transitions
			for (int i = 0; i < jStates.size(); i++) {
				// get current gene and first state of gene
				jSequence = (ArrayList) jStates.get(i);
				jState = (EmissionState) jSequence.get(0);

				// find the exo. nucl. prob. arrays based on family of Gene
				anno = jState.getAnnotation();
				si = (StateInfo) anno.getProperty(null);
				familyIndex = findFamilyIndex(si.geneName);

				// exo start array
				mean_start = probHolder.J_start_exo_mean_Fields[familyIndex - 1];
				stdDev_start = probHolder.J_start_exo_stdDev_Fields[familyIndex - 1];
				J_start_exo_probs = probHolder.getExoProbArray(mean_start,
						stdDev_start, MIN_PROB_LIMIT, false);

				// X12 to JX1 transition
				vdj.createTransition(X12_states[i], jState);

				// X11b to JXn transition (exo nucl. transitions)
				for (int j = 1; j < J_start_exo_probs.length; j++) {
					jState = (EmissionState) jSequence.get(j);
					vdj.createTransition(X11b_states[i], jState);
				}
			}

			// *** create J to J transitions *** //

			for (int i = 0; i < jStates.size(); i++) {
				jSequence = (ArrayList) jStates.get(i);

				int j;
				for (j = 0; j < (jSequence.size() - 1); j++) // for every J
															 // state but last
				{
					currJ = (EmissionState) jSequence.get(j);
					nextJ = (EmissionState) jSequence.get(j + 1);

					// J to next J transition
					vdj.createTransition(currJ, nextJ);

					// if no there was no C-region, j-gene may have been cut-off
					// (end part of JGene is not here),
					// so make it possible to have exo nucl. act. in end of gene
					// to allow for this
					if (!removed_C_Region) {
						// J to magical transition (possible to cut off most of
						// J-Gene)
						vdj.createTransition(currJ, vdj.magicalState());
					}

				}//--for(j)

				if (removed_C_Region) // if there was a C-region, make sure
									  // bridging nucl. at end of J-Gene can be
									  // cut-off
				{
					// make transition from second last J-State to magical, so
					// it can be cut-off
					vdj.createTransition(currJ, vdj.magicalState());
				}

				// make the JXlast to Magic (no cut-off in J-Gene (the whole
				// JGene is present) if this transition is chosen)
				vdj.createTransition(nextJ, vdj.magicalState());

			}//--for(i)

		} catch (Exception e) {
			throw new Error("createJTransitions(): " + e.getMessage());
		}//catch
	}// --createJTransitions()

	// ************************************************************** //

	/**
	 * set the transition (and emission) probabililites created in
	 * createJTransitions
	 */
	private void setJTransitions(SimpleMarkovModel vdj, DotState[] X11b_states,
			DotState[] X11a_states, DotState[] X12_states, ArrayList jStates,
			ProbabilityHolder probHolder, boolean removed_C_Region) //,
																	// EmissionState
																	// JC_bridge_state)
	{
		try {
			// find the exo. nucl. probability based on the family of gene
			double[] J_start_exo_probs = null;

			Annotation anno;
			StateInfo si;
			int familyIndex;
			double mean_start, stdDev_start;

			// local vars
			ArrayList jSequence; // one gene sequence
			EmissionState jState;
			EmissionState currJ = null;
			EmissionState nextJ = null;
			EmissionState lastJ = null;
			Distribution dist;

			// set the X12's to JX1's, X11b's to X11a's and X11b's to JXn+1
			// transition, plus Jx to Jx+1 probs

			double exo_prob; // use to calculate the remainding transition prob
							 // from currJ to nextJ
			// when there is exo activity involved
			double remainder; // the remainder used for JX to JX+1 trans. probs.

			int jSize; // size of dSequence
			int jGeneCount = jStates.size(); // total number of J genes

			for (int i = 0; i < jStates.size(); i++) {
				jSequence = (ArrayList) jStates.get(i);
				jSize = jSequence.size();
				jState = (EmissionState) jSequence.get(0);

				// find the exo. nucl. prob. based on family of Gene

				anno = jState.getAnnotation();
				si = (StateInfo) anno.getProperty(null);
				familyIndex = findFamilyIndex(si.geneName);

				// exo start
				mean_start = probHolder.J_start_exo_mean_Fields[familyIndex - 1];
				stdDev_start = probHolder.J_start_exo_stdDev_Fields[familyIndex - 1];
				J_start_exo_probs = probHolder.getExoProbArray(mean_start,
						stdDev_start, MIN_PROB_LIMIT, false);

				double J_start_no_exo_prob = J_start_exo_probs[0]; // the
																   // probability
																   // of zero
																   // removals

				//System.out.println("no exo prob = " + J_start_no_exo_prob);

				double X12itoJX1sProb = 1.0; // always 1.0, one to one
											 // relationship
				double X11biToX11aiProb = J_start_no_exo_prob; // 1 to 1
															   // relationship

				// set X12's to JX1's transition
				dist = vdj.getWeights(X12_states[i]);
				dist.setWeight(jState, X12itoJX1sProb); // one to many
														// relationship prob

				// set X11b's to X11a's transition
				dist = vdj.getWeights(X11b_states[i]);
				dist.setWeight(X11a_states[i], X11biToX11aiProb); // one to one
																  // relationship

				// set JXi to JXi+1
				dist = vdj.getWeights(jState);
				dist.setWeight((EmissionState) jSequence.get(1), 1.0);

				int j;
				for (j = 1; j < (jSequence.size() - 2); j++) // for every J
															 // state but first,
															 // second last and
															 // last
				{
					currJ = (EmissionState) jSequence.get(j);
					nextJ = (EmissionState) jSequence.get(j + 1);
					lastJ = (EmissionState) jSequence.get(j + 2); // this is the
																  // last J
																  // state only
																  // when the
																  // for loop
																  // has
																  // finished

					// test for exo nuclease activity at start
					if (j < J_start_exo_probs.length) { // we have exo
														// actitivity at the
														// start
						exo_prob = J_start_exo_probs[j];
						dist = vdj.getWeights(X11b_states[i]);
						dist.setWeight(currJ, exo_prob); // one to one
														 // relationship
					}
					//else // no exo actitivity

					remainder = 1.0; // the complete probability from a J-state

					if (!removed_C_Region) {
						// set the curr J to magical state transition
						double JtoMagicalProb = 0.02; // constant (not based on
													  // data)
						dist = vdj.getWeights(currJ);
						dist.setWeight(vdj.magicalState(), JtoMagicalProb);

						remainder -= JtoMagicalProb;
					}

					// set J to next J transition
					dist = vdj.getWeights(currJ);
					dist.setWeight(nextJ, remainder);

				}//--for states in curr J-Gene

				double secondLastJtoMagicalProb;
				double secondLastJtoLastJProb;

				if (removed_C_Region) {
					// set the JXlast to Magic, and JXSecondlast to Magic
					// transition
					// therefore possible to cut off last nucl. in J gene

					secondLastJtoMagicalProb = 0.5;
					secondLastJtoLastJProb = 1.0 - secondLastJtoMagicalProb; // dependant
				} else {
					// set the JXsecondLast to last J transition
					secondLastJtoMagicalProb = 0.02;
					secondLastJtoLastJProb = 1.0 - secondLastJtoMagicalProb; // dependant
				}

				// set the JXsecondLast to magical and JXsecondLast to last J
				// transition
				dist = vdj.getWeights(nextJ);
				dist.setWeight(vdj.magicalState(), secondLastJtoMagicalProb);
				dist.setWeight(lastJ, secondLastJtoLastJProb);

				// set the JXlast to Magical
				double JXlastToMagicalProb = 1.0; // always 1.0
				dist = vdj.getWeights(lastJ);
				dist.setWeight(vdj.magicalState(), JXlastToMagicalProb);

			}//--for every gene in arraylist

		} catch (Exception e) {
			throw new Error("setJTransitions:  " + e.getMessage());
		}//catch
	}// --setJTransitions()

	///////////////////////**********************************////////////////////////
	// **************************** Helper Methods
	// ********************************//
	///////////////////////**********************************////////////////////////

	/**
	 * calculate the sum of all elements in an array returns the sum of all
	 * elements in double array "array"
	 */
	private double getSum(double[] array) {
		double sum = 0;
		for (int i = 0; i < array.length; i++) {
			sum += array[i];
		}

		return sum;
	}//--getSum

	/**
	 * calculate the sum of all elements in an array returns the sum of all
	 * elements in integer array "array"
	 */
	private int getSum(int[] array) {
		int sum = 0;
		for (int i = 0; i < array.length; i++) {
			sum += array[i];
		}

		return sum;
	}//--getSum

	/**
	 * set the emission probabilities in the P-addtion states representing the
	 * possible p-region at the end of a Gene sequence
	 */
	private void set_P_end(Sequence gene, EmissionState[] pEnd) {
		int pEndLength = pEnd.length;

		String geneString = gene.seqString();

		// get the last X nucl. as a String
		String geneEnd = geneString.substring(geneString.length() - pEndLength);

		// get a nucleotide string in reversed order
		String geneEndRC = getReversedCompliment(geneEnd);

		Distribution dist;
		char nucl;
		for (int i = 0; i < pEndLength; i++) {
			// get the distribution of current p end state
			dist = pEnd[i].getDistribution();

			// get the nucleotide emitted by p end state
			nucl = geneEndRC.charAt(i);

			// set the probability of emissions
			setPEmission(dist, nucl);
			//	System.out.print(nucl);
		}//--for(i)
	}//--set_P_end()

	/**
	 * set the emission probabilities in the P-addtion states representing the
	 * possible p-region at the start of a Gene sequence
	 */
	private void set_P_start(Sequence gene, EmissionState[] pStart) {
		int pStartLength = pStart.length;

		String geneString = gene.seqString();

		// get the last X nucl. as a String
		String geneStart = geneString.substring(0, pStartLength);

		String geneStartRC = getReversedCompliment(geneStart);

		Distribution dist;
		char nucl;
		for (int i = 0; i < pStartLength; i++) {
			// get the distribution of current p start state
			dist = pStart[i].getDistribution();

			// get the nucleotide emitted by p start state
			nucl = geneStartRC.charAt(i);

			// set the probability of emissions
			setPEmission(dist, nucl);
			//	System.out.print(nucl);
		}//--for(i)
	}//--set_P_start()

	/**
	 * find and set the probability of a mutation at position "nucl_pos" into
	 * gene "seq" of type "gene_type" NB!, "nucl_pos" starts at pos 1 (first
	 * element at index 1)
	 */
	private void setGeneEmissionProb(Sequence seq, String seq_string,
			int nucl_pos, Distribution dist, FiniteAlphabet dna, int gene_type,
			int completeVGeneLength, int VGene_start_offset,
			double A_probability) {

		// calculate the probability of a mutation based on A.e^(-0,0024k).M

		// first find e^(-.0024k)
		// depends on type of gene

		double exp_mutation_prob;

		switch (gene_type) {
		case GENE_TYPE_V:
			exp_mutation_prob = ExponentialDecay.exponentialDecayVGene(
					seq_string, nucl_pos, VGene_start_offset);
			break;
		case GENE_TYPE_D:
			// becuase of antigen selection rules, the mutability score of the
			// DGene (actually the DGene and its VD and DJ junctions as well,
			// but that's irrelevant in this matter)
			// has to be multiplied by 1.5
			exp_mutation_prob = (double) 1.5
					* ExponentialDecay.exponentialDecayDGene(seq_string,
							nucl_pos, completeVGeneLength);
			break;
		case GENE_TYPE_J:
			exp_mutation_prob = ExponentialDecay.exponentialDecayJGene(
					seq_string, nucl_pos, completeVGeneLength);
			break;
		default:
			throw new Error("Gene type not V,D or J");

		}//--switch

		// get the mutability score (M)

		// get the pentanucleotide surrounding the current nucl. position from
		// the current gene sequence
		String pentaNucleotide = getPentaNucleotide(seq_string, nucl_pos);

		// calculate the mutability score of the current nucleotide based on the
		// pentanucleotide "pentanucleotide"
		double mutability_score = G_mutability_score
				.pentaNucleotideScore(pentaNucleotide);

		// calculate probability of mutation (A.e^(-0,0024k)M)
		double probability_of_mutation = A_probability * exp_mutation_prob
				* mutability_score;

		G_fostream.println("gene name: " + seq.getName() + "  nucl. position: "
				+ nucl_pos + "  Probability of mutation  = "
				+ probability_of_mutation);

		// calculate probability of no mutation
		double noMutationProb = 1.0 - probability_of_mutation;

		// get the trnucleotide from seq string and nucl pos
		String trinucleotide = getTriNucleotide(seq_string, nucl_pos); // the
																	   // trinucleotide
																	   // we
																	   // want
																	   // to
																	   // lookup

		// the nucleotide we mutate from ( the middle nucl. in the tri nucl.)
		char mutateFrom = seq_string.charAt((nucl_pos - 1));

		char nucl;
		double currMutationFraction; // a trinucleotide mutation share of the
									 // total mutation probability
		double relativeMutationProb; // the adjusted mutation probability of a
									 // trinucleotide

		// make the probability of each of the three possible mutation sum to
		// the total probability of a mutation occuring
		// but add a small probability 2% to every probability to avoid ZERO
		// probailities occuring

		// so first take 6% off the mutation probability (2% * 3 (all three
		// possible mutations)

		double six_percent_of_mutation_prob = probability_of_mutation * 0.06; // SIX
																			  // percent
		double two_percent_of_probability = probability_of_mutation * 0.02; // two
																			// percent
		double reduced_probability_of_mutation = probability_of_mutation
				- six_percent_of_mutation_prob;

		// for A,C,G,T
		for (int i = 0; i < G_mutationSpectrum.UNIQUE_NUCLEOTIDE_ARRAY.length; i++) {
			nucl = G_mutationSpectrum.UNIQUE_NUCLEOTIDE_ARRAY[i]; // a,c,g,t

			if (nucl == mutateFrom) // nucl. representing no mutation
			{
				// set probability of no mutation
				setNucleotideGeneEmission(nucl, noMutationProb, dist);
			}//--no mutation
			else // a mutation
			{
				// find probability of trinucleotide mutating into this
				// nucleotide
				currMutationFraction = G_mutationSpectrum.getTNProbability(
						trinucleotide, nucl);

				if (currMutationFraction == G_mutationSpectrum.VALUE_NO_MATCH) {
					throw new Error(
							"setGeneEmissionProb(): trinucleotide probability match lookup not found: "
									+ trinucleotide + " : " + nucl);
				}

				// make the probability of each of the three possible mutation
				// sum to
				// the total probability of a mutation occuring
				// but add a small probability 2% to every probability to avoid
				// ZERO probailities occuring
				relativeMutationProb = ((reduced_probability_of_mutation * currMutationFraction) + two_percent_of_probability);

				if (relativeMutationProb == 0)
					throw new Error("ZERO probability: " + relativeMutationProb);

				// make certain no mutation probability is ZERO
				// by adding a small number to every probability found
				if (relativeMutationProb == 0)
					G_fostream.println("************* mutation to nucl : "
							+ nucl + " = " + relativeMutationProb);
				else
					G_fostream.println("mutation to nucl : " + nucl + " = "
							+ relativeMutationProb);

				setNucleotideGeneEmission(nucl, relativeMutationProb, dist);
			}//--mutation
		}//--for(i)

	}//--setGeneEmissionProb()

	/**
	 * get the tri nucleotide centered on gene position "nucl_pos" from gene
	 * sequence "seq_string". at the edges, nucleotide will be replaced by
	 * NUCLEOTIDE_N (char 'n'
	 * 
	 * NB! "nucl_pos" counts from 1 (not zero)
	 */
	private String getTriNucleotide(String seq_string, int nucl_pos) {
		// the resulting trinucleotide
		String trinucleotide = null;

		// get the nucl string position
		int string_nucl_pos = nucl_pos - 1;

		if (string_nucl_pos == 0) // if first trinucleotide in sequence string
		{
			trinucleotide = G_mutationSpectrum.NUCLEOTIDE_N
					+ seq_string.substring(string_nucl_pos,
							(string_nucl_pos + 2));
		} else if (string_nucl_pos == (seq_string.length() - 1)) // if last
																 // trinucleotide
																 // in sequence
																 // string
		{
			trinucleotide = seq_string.substring((string_nucl_pos - 1),
					(string_nucl_pos + 1))
					+ G_mutationSpectrum.NUCLEOTIDE_N;
		} else // any trinucleotide but first or last in sequence
		{
			trinucleotide = seq_string.substring((string_nucl_pos - 1),
					(string_nucl_pos + 2));
		}

		// return result
		return trinucleotide;

	}//--getTrinucleotide()

	/**
	 * get the penta nucleotide centered on gene position "nucl_pos" from gene
	 * sequence "seq_string". At the edges (first and second / last and second
	 * last nucl. pos.), nucleotides will be replaced by U (char 'u')
	 * 
	 * NB! "nucl_pos" has first element at pos 1
	 */
	private String getPentaNucleotide(String seq_string, int nucl_pos) {
		// get the nucl string position
		int string_nucl_pos = nucl_pos - 1;
		int max_string_pos = seq_string.length() - 1;

		// resulting pentanucleotide
		String pentaNucleotide = null;

		// test the basket cases (replacement in beginning)
		if (string_nucl_pos == 0) { // first two nucl. are missing
			pentaNucleotide = "" + U_NUCLEOTIDE + U_NUCLEOTIDE
					+ seq_string.substring(0, 3);
		} else if (string_nucl_pos == 1) { // first nucl. is missing
			pentaNucleotide = "" + U_NUCLEOTIDE + seq_string.substring(0, 4);
		} else if (string_nucl_pos == max_string_pos) { // last two nucl. are
														// missing
			pentaNucleotide = ""
					+ seq_string.substring((max_string_pos - 2),
							(max_string_pos + 1)) + U_NUCLEOTIDE + U_NUCLEOTIDE;
		} else if (string_nucl_pos == (max_string_pos - 1)) { // last nucl. is
															  // missing
			pentaNucleotide = ""
					+ seq_string.substring((max_string_pos - 3),
							(max_string_pos + 1)) + U_NUCLEOTIDE;
		} else { // all nucl. are present
			pentaNucleotide = seq_string.substring((string_nucl_pos - 2),
					(string_nucl_pos + 3));
		}

		// return the pentanucleotide we found/created
		return pentaNucleotide;

	}//--getPentaNucleotide()

	//////////////////////////////////////////////////////////////////

	/**
	 * sets the emission probability of distribution "dist" to 100 percent
	 * chance of emitting nucleotide "nucl" and Zero percent probability of
	 * emitting any other nucleotides
	 */
	private void setNucleotideGeneEmission(char nucl, double probability,
			Distribution dist) {
		try {
			switch (nucl) {
			case 'a':
				dist.setWeight(DNATools.a(), probability);
				break;
			case 'c':
				dist.setWeight(DNATools.c(), probability);
				break;
			case 'g':
				dist.setWeight(DNATools.g(), probability);
				break;
			case 't':
				dist.setWeight(DNATools.t(), probability);
				break;
			}
		} catch (Exception ille) {
			throw new Error("setNucleotideGeneEmission():  "
					+ ille.getMessage());
		}
	}//--setNucleotideGeneEmission()

	/**
	 * set the emission a,c,g,t emission probability of distribution "pDist",
	 * based on the nucleotide "nucl". Basically, a P-state can only emit one
	 * symbol, so the nucleotide "nucl", has 100 percent probability, and the
	 * other nucleotides has a Zero percent chance of being emit from this state
	 */
	private void setPEmission(Distribution pDist, char nucl) {
		try {

			// set the emission probability for this P state
			if (nucl == 'a')
				pDist.setWeight(DNATools.a(), 1.0);
			else
				pDist.setWeight(DNATools.a(), 0.0);
			if (nucl == 'c')
				pDist.setWeight(DNATools.c(), 1.0);
			else
				pDist.setWeight(DNATools.c(), 0.0);
			if (nucl == 'g')
				pDist.setWeight(DNATools.g(), 1.0);
			else
				pDist.setWeight(DNATools.g(), 0.0);
			if (nucl == 't')
				pDist.setWeight(DNATools.t(), 1.0);
			else
				pDist.setWeight(DNATools.t(), 0.0);

			/*
			 * // set the emission probability for this P state
			 * pDist.setWeight(DNATools.a(), 0.25);
			 * pDist.setWeight(DNATools.c(), 0.25);
			 * pDist.setWeight(DNATools.g(), 0.25);
			 * pDist.setWeight(DNATools.t(), 0.25);
			 */
		} catch (Exception e) {
			throw new Error("setPEmission(): Can't set State Emission");
		}
	}

	/**
	 * get the reversed result of nucleotide string "seq" returns the reverse
	 * result of String "seq"
	 */
	private String getReversedCompliment(String seq) {
		char orig;
		char compliment;
		String result = "";

		// for every nucleotide (character) in nucleotide string (character
		// String) "seq"
		for (int i = 0; i < seq.length(); i++) {
			orig = seq.charAt(i);
			switch (orig) {
			case 'g':
				compliment = 'c';
				break;
			case 'c':
				compliment = 'g';
				break;
			case 'a':
				compliment = 't';
				break;
			case 't':
				compliment = 'a';
				break;
			default:
				throw new Error(
						"getReversedCompliment(): not a,c,g or t, but:  "
								+ orig);
			}
			// reverse sequence by adding to the front
			result = compliment + result;
		}

		return result; // the reversed string

	}//--getReversedCompliment()

}//class

