package immu;

/**
 * this class used to find Amino Acid Seqs.
 * Also replace old class to find STOP Codon
 */

public class CodonTranslate {

   //main for testing ONLY
	public static void main(String[] args) throws Exception {
		
		CodonTranslate myTranslate = new CodonTranslate();
        String mySeq, mySeq1, mySeq2;
        mySeq="caggtgcagctacagcagtggggcgcaggactgttgaagccttcggagaccctgcccctcacctgcggtgtctatggtgggtccttcactggtgacttctggacctggatccgccagcccccagggaagggactggagtggattggggaaatctatcaaagtggaagcaccaactacaacccgtccctcaagagtcgagtcaccatatcaatagccacgtccaagaaccaattctccctgaggctgaattctttgaccgccgcggacacggccaaatatttctgtgcgagaggcctctcgaatactgcaggtcgtcggggcccacccgctaaggctatggacgtctggggccaagggaccacggtcaccgtctcctca";
        mySeq1="ggtgcagctacagcagtggggcgcaggactgttgaagccttcggagaccctgcccctcacctgcggtgtctatggtgggtccttcactggtgacttctggacctggatccgccagcccccagggaagggactggagtggattggggaaatctatcaaagtggaagcaccaactacaacccgtccctcaagagtcgagtcaccatatcaatagccacgtccaagaaccaattctccctgaggctgaattctttgaccgccgcggacacggccaaatatttctgtgcgagaggcctctcgaatactgcaggtcgtcggggcccacccgctaaggctatggacgtctggggccaagggaccacggtcaccgtctcctca";
        mySeq2="nnnnnggtgcagctacagcagtggggcgcaggactgttgaagccttcggagaccctgcccctcacctgcggtgtctatggtgggtccttcactggtgacttctggacctggatccgccagcccccagggaagggactggagtggattggggaaatctatcaaagtggaagcaccaactacaacccgtccctcaagagtcgagtcaccatatcaatagccacgtccaagaaccaattctccctgaggctgaattctttgaccgccgcggacacggccaaatatttctgtgcgagaggcctctcgaatactgcaggtcgtcggggcccacccgctaaggctatggacgtctggggccaagggaccacggtcaccgtctcctca";
        
        //Seq1  Vstart=2,  3-2    Seq2 Vstart=-3  3-(-3)        
        
        System.out.println(myTranslate.findCodons(mySeq, 0));
        System.out.println(" "+myTranslate.findCodons(mySeq1, 1));     
        System.out.println(" "+myTranslate.findCodons(mySeq2, 6));
	}	
	
	

	public String findCodons(String sequenceString, int codonStartPos)
	{	
		String currCodon= "";
		String AminoAcid="";
	
		for(int i=codonStartPos; i<= sequenceString.length() - 3; i += 3)
		{
			currCodon = sequenceString.substring(i, i+3);	
			String curAA=TranslateAminoAcid(currCodon);
			AminoAcid+=curAA;
			if (curAA.equals("-"))  {
				AminoAcid += "STOP";
				break;
			}	
		}
		
	    return AminoAcid;
	}//findCodons()
	
	
	private String TranslateAminoAcid(String codon) {
		codon = codon.toUpperCase();
		String AA="";
		
		char ntOne, ntTwo, ntThree;
		ntOne = codon.charAt(0);
		ntTwo = codon.charAt(1);
		ntThree = codon.charAt(2);
		
		if (ntOne=='T'){
			if(ntTwo == 'T')
			{
				if(ntThree == 'T')      AA= "F";
				if(ntThree == 'C')      AA= "F";
				if(ntThree == 'A')      AA= "L";
				if(ntThree == 'G')      AA= "L";
			}
			
			if(ntTwo == 'C')
			{  
				AA= "S";
			}
			
			if(ntTwo == 'A')
			{
				if(ntThree == 'T')      AA= "Y";
				if(ntThree == 'C')      AA= "Y";
				if(ntThree == 'A')      AA= "-";
				if(ntThree == 'G')      AA= "-";
			}
			
			if(ntTwo == 'G')
			{
				if(ntThree == 'T')      AA= "C";
				if(ntThree == 'C')      AA= "C";
				if(ntThree == 'A')      AA= "-";
				if(ntThree == 'G')      AA= "W";
			}
		}
		
		if (ntOne=='C'){
			if(ntTwo == 'T')
			{
				AA= "L";
			}
			
			if(ntTwo == 'C')
			{  
				AA= "P";
			}
			
			if(ntTwo == 'A')
			{
				if(ntThree == 'T')      AA= "H";
				if(ntThree == 'C')      AA= "H";
				if(ntThree == 'A')      AA= "Q";
				if(ntThree == 'G')      AA= "Q";
			}
			
			if(ntTwo == 'G')
			{  
				AA= "R";
			}
				
		}
		
		if (ntOne=='A'){
			if(ntTwo == 'T')
			{
				if(ntThree == 'T')      AA= "I";
				if(ntThree == 'C')      AA= "I";
				if(ntThree == 'A')      AA= "I";
				if(ntThree == 'G')      AA= "M";
			}
			
			if(ntTwo == 'C')
			{  
				AA= "T";
			}
			
			if(ntTwo == 'A')
			{
				if(ntThree == 'T')      AA= "N";
				if(ntThree == 'C')      AA= "N";
				if(ntThree == 'A')      AA= "K";
				if(ntThree == 'G')      AA= "K";
			}
			
			if(ntTwo == 'G')
			{
				if(ntThree == 'T')      AA= "S";
				if(ntThree == 'C')      AA= "S";
				if(ntThree == 'A')      AA= "R";
				if(ntThree == 'G')      AA= "R";
			}
		}
		
		if (ntOne=='G'){
			if(ntTwo == 'T')
			{   
				AA= "V";
			}
			
			if(ntTwo == 'C')
			{  
				AA= "A";
			}
			
			if(ntTwo == 'A')
			{
				if(ntThree == 'T')      AA= "D";
				if(ntThree == 'C')      AA= "D";
				if(ntThree == 'A')      AA= "E";
				if(ntThree == 'G')      AA= "E";
			}
			
			if(ntTwo == 'G')
			{
				AA= "G";
			}
		}
		
		return AA;
	}   //TranslateAminoAcid
			
}
