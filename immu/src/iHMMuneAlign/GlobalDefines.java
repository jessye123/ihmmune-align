package iHMMuneAlign;

/**
 * @author harris
 * 
 * TODO
 *  
 */
interface GlobalDefines {
	// defines
	final char V_STATE_TOKEN = 'V';

	final char D_STATE_TOKEN = 'D';

	final char J_STATE_TOKEN = 'J';

	final char DOT_STATE_TOKEN = 'X';

	final char VD_N_TOKEN = '1';

	final char DJ_N_TOKEN = '2';

	final char V_END_P_TOKEN = '3';

	final char D_START_P_TOKEN = '4';

	final char D_END_P_TOKEN = '5';

	final char J_START_P_TOKEN = '6';

	final char JC_BRIDGE_TOKEN = 'B';

	final char N_TOKEN = 'N';

	final char P_TOKEN = 'P';

	final char U_NUCLEOTIDE = 'u';

	final char J_END_C_STATE_TOKEN = 'C';

	final char P_STATE_TOKEN = 'P';

	// the result of comparing the UMS sequence with the pre emission of a state
	// can be a
	// match, a mismatch or NA
	final char MATCH_TOKEN = '.';

	final char MISMATCH_TOKEN = '|';

	final char NA_TOKEN = '-';

	final char ANY_PRE_EMISSION_TOKEN = '?';

	final char NO_PRE_EMISSION_TOKEN = ' ';

	final char GAP_TOKEN = ' ';

	final char HTML_WRITER_NONE = 0;

	final char HTML_WRITER_EXCEL = 1;

	final char HTML_WRITER_WWW = 2;

	// alignment types, specifying what type of alignment is being performed
	final byte NO_ALIGNMENT = 0;

	final byte SINGLE_ALIGNMENT = 1;

	final byte MULTIPLE_ALIGNMENT = 2;

	final byte SINGLE_FILE_MULTIPLE_ALIGNMENT = 4;

	// standard DGene mutation/match lengths
	
	final int NO_SELECTED_DGENE_ACCEPTANCE = -1;
	
	// In order to be a proper EIGHT_MER_DGENE_ACCEPTANCE Alignment, we need minimum:
	// 8 consecutive matches with no mismatches
	// or 10 consecutive matches with no more than one mismatch
	// or 12 or more consecutive matches with up to two mismatches	
	final int EIGHT_MER_DGENE_ACCEPTANCE = 1;
	
	// In order to be a FIVE_MER_DGENE_ACCEPTANCE Alignment, we need minimum:
	// 5 consecutive matches with no mismatches	
	final int FIVE_MER_CONSECUTIVE_DGENE_ACCEPTANCE = 2;
	
	
	// old version
	char N3_STATE_TOKEN = '3';

	char N2_STATE_TOKEN = '2';

	char N1_STATE_TOKEN = '1';
	
	
}