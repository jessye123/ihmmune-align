package immu;
//imports

public class IHMMuneResult{
	//array list to hold the iHMMune-align result for each sequence from a specified file
	//ArrayList iHMMuneAlignResults = new ArrayList();
	
	//vars to hold the features of an iHMMune-align partitioning result
	private String seqID;
	private String IGHVName;
	private String IGHDName;
	private String IGHJName;
	private String IGHVSeq;
	private String N1Seq;
	private String IGHDSeq;
	private String N2Seq;
	private String IGHJSeq;
	private String VMut;
	private String DMut;
	private String JMut;
	
	//brief iHMMune-align representation with V/D/J names, seqs and mutations only
	public IHMMuneResult (String seqID, String IGHVName, String IGHDName, String IGHJName, String IGHVSeq, String N1Seq, String IGHDSeq, String N2Seq, String IGHJSeq,
							   String VMut, String DMut, String JMut) 
	{
		this.seqID = seqID;
		this.IGHVName = IGHVName;
		this.IGHDName = IGHDName;
		this.IGHJName = IGHJName;
		this.IGHVSeq = IGHVSeq;
		this.N1Seq = N1Seq;
		this.IGHDSeq = IGHDSeq;
		this.N2Seq = N2Seq;
		this.IGHJSeq = IGHJSeq;
		this.VMut = VMut;
		this.DMut = DMut;
		this.JMut = JMut;
	}
	
	//add a function to get the seqID
	public String getSeqID()
	{
		return seqID;
	}
	
	//add a function to get the IGHVName
	public String getIGHVName()
	{
		return IGHVName;
	}
	
	//get the IGHD gene name
	public String getIGHDName()
	{
		return IGHDName;
	}
	
	//get the IGHJ gene name
	public String getIGHJName()
	{
		return IGHJName;
	}
	
	//get the V region sequence
	public String getVSeq()
	{
		//char[] IGHVSeqArr = IGHVSeq.toCharArray();
		return IGHVSeq;
	}
	
	//get the N1 sequence
	public String getN1Seq() 
	{
		//char[] N1SeqArr = N1Seq.toCharArray();
		return N1Seq;
	}
	
	//get the D region seq
	public String getDSeq()
	{
		//char[] DSeqArr = IGHDSeq.toCharArray();
		return IGHDSeq;
	}
	
	//get the N2 region seq
	public String getN2Seq()
	{
		//char[] N2SeqArr = N2Seq.toCharArray();
		return N2Seq;
	}
	
	//J region sequence
	public String getJSeq() 
	{
		//char[] JSeqArr = IGHJSeq.toCharArray();
		return IGHJSeq;
	}
	
	//V mutation number
	public int getVMut()
	{
		int VMutInt = Integer.parseInt(VMut);
		//return VMutInt;
		return VMutInt;
	}
	
	//D mutation number
	public int getDMut()
	{
		int DMutInt = Integer.parseInt(DMut);
		return DMutInt;
	}
	
	//J mut no
	public int getJMut()
	{
		int JMutInt = Integer.parseInt(JMut);
		return JMutInt;
	}
	

}