package au.edu.uq.imb.memesuite.db;

public class Category {
  private long id;
  private String name;
  private long cnt;

  public Category(long id, String name, long cnt) {
    this.id = id;
    this.name = name;
    this.cnt = cnt;
  }

  public long getID() {
    return this.id;
  }

  public String getName() {
    return this.name;
  }

  public long getCnt() {
    return this.cnt;
  }
}
