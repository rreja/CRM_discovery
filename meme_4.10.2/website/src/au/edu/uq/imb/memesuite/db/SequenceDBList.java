package au.edu.uq.imb.memesuite.db;

import au.edu.uq.imb.memesuite.data.Alphabet;
import org.sqlite.SQLiteDataSource;

import java.io.File;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static au.edu.uq.imb.memesuite.db.SQL.*;

/**
 * Access information on sequence databases from the sqlite database.
 */
public class SequenceDBList extends DBList {
  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web");

  public SequenceDBList(File db) throws ClassNotFoundException {
    super(db, false);
  }

  public SequenceDBList(SQLiteDataSource ds) {
    super(ds);
  }

  @Override
  protected PreparedStatement prepareCategoriesQuery(Connection conn, boolean shortOnly, EnumSet<Alphabet> allowedAlphabets) throws SQLException {
    PreparedStatement ps = conn.prepareStatement(SQL.SELECT_SEQUENCE_CATEGORIES);
    ps.setInt(1, SQL.enumsToInt(allowedAlphabets));
    ps.setBoolean(2, shortOnly);
    return ps;
  }

  @Override
  protected PreparedStatement prepareListingsQuery(Connection conn, long categoryId, boolean shortOnly, EnumSet<Alphabet> allowedAlphabets) throws SQLException {
    PreparedStatement ps = conn.prepareStatement(SELECT_SEQUENCE_LISTINGS_OF_CATEGORY);
    ps.setLong(1, categoryId);
    ps.setInt(2, SQL.enumsToInt(allowedAlphabets));
    ps.setBoolean(3, shortOnly);
    return ps;
  }

  private SequenceDB loadSequenceFile(ResultSet rset, String listingName, String listingDescription) throws SQLException {
    long id = rset.getLong(1);
    Alphabet alphabet = Alphabet.fromInt(rset.getInt(2));
    long edition = rset.getLong(3);
    String version = rset.getString(4);
    String description = rset.getString(5);
    String fileSeq = rset.getString(6);
    String fileBg = rset.getString(7);
    long sequenceCount = rset.getLong(8);
    long totalLen = rset.getLong(9);
    long minLen = rset.getLong(10);
    long maxLen = rset.getLong(11);
    double avgLen = rset.getDouble(12);
    double stdDLen = rset.getDouble(13);
    return new SequenceDB(listingName, listingDescription, id, alphabet,
        edition, version, description, fileSeq, fileBg, sequenceCount,
        totalLen, minLen, maxLen, avgLen, stdDLen);
  }

  public List<SequenceVersion> getVersions(long listingId) throws SQLException {
    return getVersions(listingId, false, EnumSet.allOf(Alphabet.class));
  }

  public List<SequenceVersion> getVersions(long listingId, boolean shortOnly, EnumSet<Alphabet> allowedAlphabets) throws SQLException {
    List<SequenceVersion> items = new ArrayList<SequenceVersion>();
    Connection conn = null;
    PreparedStatement stmt = null;
    ResultSet rset = null;
    try {
      conn = ds.getConnection();
      // lookup up the listing
      stmt = conn.prepareStatement(SELECT_LISTING);
      stmt.setLong(1, listingId);
      rset = stmt.executeQuery();
      if (!rset.next()) throw new SQLException("No listing found for identifier " + listingId);
      final String listingName = rset.getString(1);
      final String listingDescription = rset.getString(2);
      rset.close(); rset = null;
      stmt.close(); stmt = null;
      // now lookup the sequence files and cluster them by the version
      stmt = conn.prepareStatement(SELECT_SEQUENCE_FILES_OF_LISTING);
      stmt.setLong(1, listingId);
      stmt.setInt(2, SQL.enumsToInt(allowedAlphabets));
      stmt.setBoolean(3, shortOnly);
      rset =  stmt.executeQuery();
      List<SequenceDB> files = new ArrayList<SequenceDB>();
      while (rset.next()) {
        SequenceDB file = loadSequenceFile(rset, listingName, listingDescription);
        if (!files.isEmpty() && files.get(0).getEdition() != file.getEdition()) {
          items.add(new SequenceVersion(files.get(0).getEdition(), files.get(0).getVersion(), files));
          files = new ArrayList<SequenceDB>();
        }
        files.add(file);
      }
      if (!files.isEmpty()) {
        items.add(new SequenceVersion(files.get(0).getEdition(), files.get(0).getVersion(), files));
      }
      rset.close(); rset = null;
      stmt.close(); stmt = null;
      conn.close(); conn = null;
    } finally {
      if (rset != null) {
        try {
          rset.close();
        } catch (SQLException e) { /* ignore */ }
      }
      if (stmt != null) {
        try {
          stmt.close();
        } catch (SQLException e) { /* ignore */ }
      }
      if (conn != null) {
        try {
          conn.close();
        } catch (SQLException e) { /* ignore */ }
      }
    }
    return items;
  }

  public SequenceDB getSequenceFile(long listingId, long edition, EnumSet<Alphabet> alphabets) throws SQLException {
    Connection conn = null;
    PreparedStatement pstmt = null;
    ResultSet rset = null;
    SequenceDB sequenceDB = null;
    try {
      conn = ds.getConnection();
      // get the name and description
      pstmt = conn.prepareStatement(SQL.SELECT_LISTING);
      pstmt.setLong(1, listingId);
      rset = pstmt.executeQuery();
      if (!rset.next()) throw new SQLException("No listing by that ID");
      String listingName = rset.getString(1);
      String listingDescription = rset.getString(2);
      rset.close(); rset = null;
      pstmt.close(); pstmt = null;
      pstmt = conn.prepareStatement(SELECT_SEQUENCE_FILE_OF_LISTING);
      pstmt.setLong(1, listingId);
      pstmt.setLong(2, edition);
      pstmt.setInt(3, SQL.enumsToInt(alphabets));
      rset = pstmt.executeQuery();
      if (!rset.next()) throw new SQLException("No sequence file by that combination of listing, edition and alphabets");
      sequenceDB = loadSequenceFile(rset, listingName, listingDescription);
      rset.close(); rset = null;
      pstmt.close(); pstmt = null;
      conn.close(); conn = null;
    } finally {
      if (rset != null) {
        try {
          rset.close();
        } catch (SQLException e) { /* ignore */ }
      }
      if (pstmt != null) {
        try {
          pstmt.close();
        } catch (SQLException e) { /* ignore */ }
      }
      if (conn != null) {
        try {
          conn.close();
        } catch (SQLException e) { /* ignore */ }
      }
    }
    return sequenceDB;
  }
}
