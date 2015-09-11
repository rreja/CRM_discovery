package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.Alphabet;
import au.edu.uq.imb.memesuite.db.*;
import au.edu.uq.imb.memesuite.servlet.util.WebUtils;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLSubGenerator;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.template.HTMLTemplateCache;

import javax.servlet.ServletContext;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.Collections;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.CACHE_KEY;
import static au.edu.uq.imb.memesuite.servlet.ConfigurationLoader.SEQUENCE_DB_KEY;

/**
 * Display the available sequence databases.
 */
public class ShowSequenceDBs extends HttpServlet {
  private ServletContext context;
  private HTMLTemplate template;
  private HTMLTemplate categoryTemplate;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web.sequencedb");

  public ShowSequenceDBs() { }

  @Override
  public void init() throws ServletException {
    this.context = this.getServletContext();
    HTMLTemplateCache cache = (HTMLTemplateCache)context.getAttribute(CACHE_KEY);
    template = cache.loadAndCache("/WEB-INF/templates/show_sequence_dbs.tmpl");
    categoryTemplate = template.getSubtemplate("category");
  }

  @Override
  public void doGet(HttpServletRequest request, HttpServletResponse response)
      throws IOException, ServletException {
    if (request.getParameter("category") != null) {
      outputXmlListingsOfCategory(response, getId(request, "category"),
          isShortOnly(request), allowedAlphabets(request));
    } else if (request.getParameter("listing") != null) {
      outputXmlVersionsOfListing(response, getId(request, "listing"),
          isShortOnly(request), allowedAlphabets(request));
    } else {
      display(response);
    }
  }

  @Override
  public void doPost(HttpServletRequest request, HttpServletResponse response)
      throws IOException, ServletException {
    display(response);
  }

  private void display(HttpServletResponse response)
      throws IOException, ServletException {
    SequenceDBList sequenceDBList = (SequenceDBList)context.getAttribute(SEQUENCE_DB_KEY);
    response.setContentType("text/html; charset=UTF-8");
    HTMLSub out = template.toSub();
    if (sequenceDBList != null) {
      try {
        out.set("category", new AllCategory(sequenceDBList));
      } catch (SQLException e) {
        logger.log(Level.SEVERE, "Error loading sequence categories", e);
        out.empty("category");
      }
    } else {
      out.empty("category");
    }
    out.output(response.getWriter());
  }

  private long getId(HttpServletRequest request, String name) throws ServletException {
    String value = request.getParameter(name);
    if (value == null) {
      throw new ServletException("Parameter '" + name + "' was not set.");
    }
    long id;
    try {
      id = Long.parseLong(value, 10);
    } catch (NumberFormatException e) {
      throw new ServletException("Parameter '" + name + "' is not a integer value.", e);
    }
    return id;
  }

  private boolean isShortOnly(HttpServletRequest request) throws ServletException{
    String value = request.getParameter("short");
    if (value == null) return false;
    boolean shortOnly;
    try {
      shortOnly = Integer.parseInt(value, 2) != 0;
    } catch (NumberFormatException e) {
      throw new ServletException("Parameter 'short' is not a binary value.", e);
    }
    return shortOnly;
  }

  private EnumSet<Alphabet> allowedAlphabets(HttpServletRequest request) throws ServletException {
    String value = request.getParameter("alphabets");
    if (value == null) return EnumSet.allOf(Alphabet.class);
    EnumSet<Alphabet> allowedAlphabets;
    try {
      allowedAlphabets = SQL.intToEnums(Alphabet.class, Integer.parseInt(value, 10));
    } catch (NumberFormatException e) {
      throw new ServletException("Parameter 'alphabets' is not a bitset.", e);
    } catch (IllegalArgumentException e) {
      throw new ServletException("Parameter 'alphabets' is not a bitset", e);
    }
    return allowedAlphabets;
  }

  private void outputXmlListingsOfCategory(HttpServletResponse response,
      long categoryId, boolean shortOnly, EnumSet<Alphabet> allowedAlphabets)
      throws IOException, ServletException {
    SequenceDBList sequenceDBList = (SequenceDBList)context.getAttribute(SEQUENCE_DB_KEY);
    response.setContentType("application/xml; charset=UTF-8");
    // turn off caching
//    response.setHeader("Cache-Control", "no-cache, no-store, must-revalidate"); // HTTP 1.1.
//    response.setHeader("Pragma", "no-cache"); // HTTP 1.0.
//    response.setDateHeader("Expires", 0); // Proxies.
    PrintWriter out = null;
    try {
      Iterator<Listing> listingView = sequenceDBList.getListings(categoryId, shortOnly, allowedAlphabets).iterator();
      out = response.getWriter();
      out.println("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
      out.println("<listings category=\"" + categoryId + "\" short=\"" +
          (shortOnly ? "1" : "0") +"\" alphabets=\"" +
          SQL.enumsToInt(allowedAlphabets) + "\">");
      while (listingView.hasNext()) {
        Listing listing = listingView.next();
        out.println("<l i=\"" + listing.getID() + "\" n=\"" + WebUtils.escapeForXML(listing.getName()) +
            "\" a=\"" + SQL.enumsToInt(listing.getAlphabets()) + "\"/>");
      }
      out.println("</listings>");
    } catch (SQLException e) {
      throw new ServletException(e);
    } finally {
      if (out != null) out.close();
    }
  }

  private void outputXmlVersionsOfListing(HttpServletResponse response,
      long listingId, boolean shortOnly, EnumSet<Alphabet> allowedAlphabets)
      throws IOException, ServletException {
    SequenceDBList sequenceDBList = (SequenceDBList)context.getAttribute(SEQUENCE_DB_KEY);
    response.setContentType("application/xml; charset=UTF-8");
    PrintWriter out = null;
    try {
      List<SequenceVersion> versions = sequenceDBList.getVersions(listingId, shortOnly, allowedAlphabets);
      String listingName = "";
      String listingDescription = "";
      if (!versions.isEmpty()) {
        SequenceVersion ver = versions.get(0);
        listingName = ver.getListingName();
        listingDescription = ver.getListingDescription();
      }
      out = response.getWriter();
      out.println("<versions listing=\"" + listingId + "\" name=\"" + WebUtils.escapeForXML(listingName) +
          "\" description=\"" + WebUtils.escapeForXML(listingDescription) + "\" short=\"" +
          (shortOnly ? "1" : "0") + "\" alphabets=\"" + SQL.enumsToInt(allowedAlphabets) + "\">");
      for (SequenceVersion version : versions) {
        out.println("<version id=\"" + version.getEdition() + "\" name=\"" +
            version.getVersion() + "\" alphabets=\"" +
            SQL.enumsToInt(version.getAlphabets()) + "\">");
        for (Alphabet alphabet : version.getAlphabets()) {
          SequenceDB file = version.getSequenceFile(alphabet);
          out.println("<file" +
              " alphabet=\"" + SQL.enumsToInt(Collections.singleton(file.guessAlphabet())) + "\"" +
              " count=\"" + file.getSequenceCount() + "\"" + " min=\"" + file.getMinLength() + "\"" +
              " max=\"" + file.getMaxLength() + "\"" + " avg=\"" + file.getAverageLength() + "\"" +
              " total=\"" + file.getTotalLength() + "\"" +
              " description=\"" + WebUtils.escapeForXML(file.getDescription()) + "\"" +
              "/>");
        }
        out.println("</version>");
      }
      out.println("</versions>");
    } catch (SQLException e) {
      throw new ServletException(e);
    } finally {
      if (out != null) out.close();
    }

  }

  private class AllCategory extends HTMLSubGenerator<Category> {
    private AllCategory(SequenceDBList sequenceDBList) throws SQLException{
      super(sequenceDBList.getCategories(false, EnumSet.allOf(Alphabet.class)));
    }

    @Override
    protected HTMLSub transform(Category item) {
      HTMLSub out = categoryTemplate.toSub();
      out.set("id", item.getID());
      out.set("name", item.getName());
      out.set("cnt", item.getCnt());
      return out;
    }
  }

}
