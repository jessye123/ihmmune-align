package iHMMuneAlign;

/**
 * @author harris
 * 
 * TODO
 *  
 */
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;

public class FastaReader {

	/**
	 * for calculating the average DGene length
	 */
	public static void main(String[] args) {
		FastaReader fastaReader = new FastaReader();
		try {
			ArrayList dGenes = fastaReader.readFile(new java.io.File("src/IGHD_Repertoire.fa"));

			int totalLength = 0; // combined length of all the DGenes
			int numberOfDGenes = dGenes.size();

			Sequence currGene; // current gene from list

			// for every gene in list
			for (int i = 0; i < numberOfDGenes; i++) {
				currGene = (Sequence) dGenes.get(i);

				// add the length of this gene to the total of all genes
				totalLength += currGene.length();

			}//--for()

			// find the average gene length
			double averageLength = (double) totalLength
					/ (double) numberOfDGenes;

			System.out.println("average D-Gene length = " + averageLength);
		} catch (Exception ex) {
			throw new Error(ex.getMessage());
		}//--catch

	}//--main()

	public ArrayList readFile(java.io.File file){
		ArrayList sequences = new ArrayList();
		// Set up sequence iterator
		
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(file));
		} catch(FileNotFoundException fnfe)
		{
			System.out.println("File not found");
			return null;
		}
		SequenceIterator stream = SeqIOTools.readFastaDNA(br);

		// Iterate over all sequences in the stream
		int counter = 0;
		Sequence seq = null;
		while (stream.hasNext()) {
			System.out.println(" begun reading sequence: " + counter);
			try {
				seq = stream.nextSequence();
				System.out.println(seq.toString() +"  "+ seq.seqString());
			} catch(BioException bioe)
			{
				System.out.println("Bio Exception when attempting to read from fasta file");
				return null;
			}
			sequences.add(seq);
			counter++;
		}//while
		return sequences;
	}// --f()
}// --class

