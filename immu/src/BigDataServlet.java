

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import javax.servlet.ServletException;
import javax.servlet.annotation.MultipartConfig;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.http.Part;

/**
 * Servlet implementation class BigDataServlet
 */
@WebServlet("/BigDataServlet")
@MultipartConfig(fileSizeThreshold=1024*1024*10, 	// 10 MB 
maxFileSize=1024*1024*50,      	// 50 MB
maxRequestSize=1024*1024*100)   	// 100 MB

public class BigDataServlet extends HttpServlet {
	private static final long serialVersionUID = 1L;
       
    /**
     * @see HttpServlet#HttpServlet()
     */
    public BigDataServlet() {
        super();
        // TODO Auto-generated constructor stub
    }

	/**
	 * @see HttpServlet#doGet(HttpServletRequest request, HttpServletResponse response)
	 */
	protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		// TODO Auto-generated method stub
		response.getWriter().append("Served at: ").append(request.getContextPath());
	}

	/**
	 * @see HttpServlet#doPost(HttpServletRequest request, HttpServletResponse response)
	 */
	protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
		//String uploadFilePath =   "/srvr/z2283813/tomcat7/webapps/immuFile";
		String uploadFilePath =   "/srvr/ihmmune/tomcat7/webapps/immuFile";
	       /*
	        File fileSaveDir = new File(uploadFilePath);
	        if (!fileSaveDir.exists()) {
	            fileSaveDir.mkdirs();
	        }  */
	        String type[] = request.getParameterValues("type");
	    	Part part= request.getPart("file");
	    	//String fileName = (part.getName()).toString();
	    	String fileName = "outputx9x9x9.csv";
	        if(fileName  != null && !fileName.isEmpty()) {     
	        part.write(uploadFilePath + File.separator + fileName);	           		
	        }
	        
	        if (type != null && type.length != 0) {
	        	response.getWriter().append("You have selected: ");
	        	for (int i = 0; i < type.length; i++) {
	        		response.getWriter().append(type[i]); 
	        		 request.setAttribute(type[i], type[i]);
	        	}
	        	}
	    	
	        
	        request.setAttribute("message", fileName + " File uploaded successfully!");
	        request.setAttribute("filename", fileName);
	       // request.setAttribute("type", type[]);	
			getServletContext().getRequestDispatcher("/GraphicalTool.jsp").forward(request, response); 
	      
		
	}

}
