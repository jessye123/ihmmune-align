package iHMMuneAlign;


/**
 * Rewrite By Jessica Ye
 * @author harris
 */

/**
 * holds information about one identified and aligned gene region
 */
public class GeneInfo {
	int start_gene_pos; // the alignment start position (counting from 1)
	int end_gene_pos;
	int gene_length; // length of the complete gene
	String gene_name; // name of the gene
	String aligned_gene_string;
	String aligned_ums_string;
	Boolean acceptedAlignment;
	
	/**
	 * constructor
	 */
	public GeneInfo(String geneName, String geneString, String umsString,
			int geneLength, int startGenePos, int endGenePos) {
		this.gene_name = geneName;
		this.aligned_gene_string = geneString;
		this.aligned_ums_string = umsString;
		this.gene_length = geneLength;
		this.start_gene_pos = startGenePos; // the position of the first aligned
											// nucl. (counting from pos 1)
		this.end_gene_pos = endGenePos; // the position of the last aligned
										// nucl. (counting from pos 1)
	}//--GeneInfo
	
	public int getMutation(){
		String Gene = aligned_gene_string.toLowerCase();
		String matching_sequence = aligned_ums_string.toLowerCase();
		int number_of_mutations_in_Gene = 0; 
		int number_of_NandX_nts_in_Gene = 0; 
		char V_nucl, MS_nucl; // current VGene and mathcing sequence Nucleotide

		for (int i = 0; i < Gene.length(); i++) {
				V_nucl = Gene.charAt(i);
				MS_nucl = matching_sequence.charAt(i);			
				if (V_nucl != MS_nucl) {				
					if (MS_nucl == 'N' || MS_nucl == 'X') {
						++number_of_NandX_nts_in_Gene; 
					} else {					
						++number_of_mutations_in_Gene; 
					}
				}
			}//--for(i)
	return number_of_mutations_in_Gene;
	}

}//--GeneInfo
