package iHMMuneAlign;

/**
 * @author harris
 *
 * TODO
 * 
 */
import java.util.*;
import org.biojava.bio.*;
import org.biojava.utils.*;

/**
 * For each state in a sequence, this class mimics an Annotation object only to
 * hold certain information about that state
 */
public class FakeAnnotation implements Annotation {
	// holds all the info on this state
	StateInfo stateInfo = null;

	public FakeAnnotation(StateInfo stateInfo) {
		this.stateInfo = stateInfo;
	}

	public static void main(String[] args) {
		StateInfo siTest = new StateInfo("abc", 'D', 10, 8, 'e', -1);
		FakeAnnotation fa = new FakeAnnotation(siTest);
		if (!fa.containsProperty(null))
			fa.setProperty(null, siTest);
		StateInfo si = (StateInfo) fa.getProperty(null);

		System.out.println("si nucloetide number = "
				+ si.getNucleotidePosition());

	}

	public Map asMap() {
		return null;
	}

	public boolean containsProperty(Object key) {
		if (stateInfo == null)
			return false;
		else
			return true;
	}

	public Object getProperty(Object key) {
		return stateInfo;
	}

	public Set keys() {
		return null;
	}

	public void removeProperty(Object key) {
		stateInfo = null;
	}

	public void setProperty(Object key, Object value) {
		stateInfo = (StateInfo) value;
	}

	public void addChangeListener(ChangeListener cl) {
	}

	public void addChangeListener(ChangeListener cl, ChangeType ct) {
	}

	public void removeChangeListener(ChangeListener cl) {
	}

	public void removeChangeListener(ChangeListener cl, ChangeType ct) {
	}

	public boolean isUnchanging(ChangeType ct) {
		return true;
	}
}

class StateInfo implements GlobalDefines {
	//class members
	String geneName;

	char geneToken;

	char preEmissionSymbol;

	int geneNumber;

	int nucleotidePosition;

	int completeGeneLength;

	public StateInfo(String geneName, char geneToken, int geneNumber,
			int nucleotidePosition, char preEmissionSymbol,
			int completeGeneLength) {
		this.geneName = geneName;
		this.geneToken = geneToken;
		this.geneNumber = geneNumber;
		this.nucleotidePosition = nucleotidePosition;
		this.preEmissionSymbol = preEmissionSymbol;
		this.completeGeneLength = completeGeneLength;
	}

	public StateInfo(String geneName, char geneToken) {
		this.geneName = geneName;
		this.geneToken = geneToken;
		// if this is a dot state
		if (geneToken == DOT_STATE_TOKEN)
			this.preEmissionSymbol = NO_PRE_EMISSION_TOKEN;
		else
			this.preEmissionSymbol = ANY_PRE_EMISSION_TOKEN;
	}

	public String getGeneName() {
		return geneName;
	}

	public char getGeneToken() {
		return geneToken;
	}

	public int getGeneNumber() {
		return geneNumber;
	}

	public int getNucleotidePosition() {
		return nucleotidePosition;
	}

	public char getPrepreEmissionSymbol() {
		return preEmissionSymbol;
	}
}

/*
 * class StateSet extends AbstractSet { public int size() { return 3; }
 * 
 * public Iterator iterator() { return null; } }
 *  
 */
