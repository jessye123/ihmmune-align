
import java.io.File;
import java.io.PrintWriter;
import java.io.IOException;
import javax.servlet.ServletException;
import javax.servlet.annotation.MultipartConfig;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.http.Part;

/**
 * Servlet implementation class Mytest
 */
@WebServlet("/MultiServlet")
@MultipartConfig(fileSizeThreshold=1024*1024*10, 	// 10 MB 
maxFileSize=1024*1024*50,      	// 50 MB
maxRequestSize=1024*1024*100)   	// 100 MB

public class MultiServlet extends HttpServlet {
	private static final long serialVersionUID = 1L;
       
    /**
     * @see HttpServlet#HttpServlet()
     */
    public MultiServlet() {
        super();
        // TODO Auto-generated constructor stub
    }

	/**
	 * @see HttpServlet#doGet(HttpServletRequest request, HttpServletResponse response)
	 */
	protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		// TODO Auto-generated method stub
		response.getWriter().append("Served at: ").append(request.getContextPath());
		(response.getWriter()).println("<h1>Hello, World!</h1>");
	}

	/**
	 * @see HttpServlet#doPost(HttpServletRequest request, HttpServletResponse response)
	 */
	protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		      
        String uploadFilePath =   "/srvr/ihmmune/tomcat7/webapps/immuFile";
       /*
        File fileSaveDir = new File(uploadFilePath);
        if (!fileSaveDir.exists()) {
            fileSaveDir.mkdirs();
        }  */
   
        String type = request.getParameter("type");
        
    	Part part= request.getPart("file");
        String fileName = (part.getName()).toString();// part.getSubmittedFileName();
        
        if(fileName  != null && !fileName.isEmpty()) {     
        part.write(uploadFilePath + File.separator + fileName);
           		
        }
    	else {
    		 String dna = request.getParameter("DNA");
    		 fileName="myfile.txt";
    		 //try(  PrintWriter out = new PrintWriter(uploadFilePath+ File.separator + fileName )  ){
    			////    out.println( dna.trim() );
    			//}
    		(response.getWriter()).println("<h1>No file! DNA String:</h1>"+dna +type);
    	}
        
        request.setAttribute("message", fileName + " File uploaded successfully!");
        request.setAttribute("filename", fileName);
        request.setAttribute("type", type);	
		getServletContext().getRequestDispatcher("/MultiSeqTool.jsp").forward(request, response); 
      
	}

}
