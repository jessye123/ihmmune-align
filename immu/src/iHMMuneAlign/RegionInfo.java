package iHMMuneAlign;

/**
 * Rewrite By Jessica Ye
 * @author harris
 *
 * holds information about one identified and aligned N or P region 
 */

public class RegionInfo {
	int start_pos; // the sequence position
	int end_pos;
	String region_string;// the nucleotide sequence making up this region

	/**
	 * constructor
	 */
	public RegionInfo(String regionString, int startPos, int endPos) {
		this.region_string = regionString;
		this.start_pos = startPos;
		this.end_pos = endPos;

	}//--RegionInfo()

}//--RegionInfo
