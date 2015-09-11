package au.edu.uq.imb.memesuite.tomcat.valves;

// This valve limits restrict the number of jobs users (as identified by IP address)
// can submitt in a time interval. Addresses provided in a 'whitelist' are
// exmpt from restrictions.

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.ServletException;
import org.apache.catalina.connector.Request;
import org.apache.catalina.connector.Response;
import org.apache.catalina.LifecycleException;
import org.apache.catalina.Valve;
import org.apache.catalina.valves.ValveBase;

// This class tracks the number of jobs submitted from an IP address
// If the number of jobs from the IP reaches the jobLimit parameter
// no more jobs be added. Jobs that are older then the duration
// parameter (in seconds) are deleted from the recored for the IP
// address. When all the jobs for an address have deleted, the record
// for that address will be deleted.
class TimeLimitedJobTable {

    private Hashtable<String, Vector<Date> > jobTable;
    private long duration; // In seconds
    private int jobLimit;

    TimeLimitedJobTable(long duration, int jobLimit) {
        this.jobTable = new Hashtable<String, Vector<Date> >();
        this.duration = duration;
        this.jobLimit = jobLimit;
    }

    // Add a job to the job table unless the user has
    // too many jobs in the time interval.
    synchronized boolean addJob(String userAddress) {
      boolean succeded = false;
      flushOldJobs();
      Vector<Date> jobs = jobTable.get(userAddress);
      if (jobs == null) {
        jobs = new Vector<Date>();
        jobTable.put(userAddress, jobs);
      }
      if (jobs.size() < jobLimit) {
        succeded = jobs.add(new Date());
      }
      return succeded;
    }

    synchronized int getNumJobs(String userAddress) {
      int numJobs = 0;
      if (jobTable.containsKey(userAddress)) {
        numJobs = jobTable.get(userAddress).size();
      }
      return numJobs;
    }

    // Calculate how long the user should wait before
    // submitting more jobs. 
    // Return 0 if the user has fewer jobs then the 
    // maximum allowed.
    synchronized long getMinWait(String userAddress) {
      long minWait = 0;
      Vector<Date> jobs = jobTable.get(userAddress);
      if (jobs != null && jobs.size() >= jobLimit) {
        Date oldestDate = jobs.firstElement();
        Date currentDate = new Date();
        long age = (currentDate.getTime() - oldestDate.getTime()) / 1000;
        if (age < duration) {
          minWait = duration - age;
        }
      }
      return minWait;
    }

    // Iterate over the table removing records that have expired.
    synchronized void flushOldJobs() {
      Date currentDate = new Date();
      // Loop over all known user addreses
      Enumeration<String> e = jobTable.keys();
      while(e.hasMoreElements()) {
        String key = e.nextElement();
        Vector<Date> jobs = jobTable.get(key);
        // Loop over all jobs for the user address
        Iterator<Date> it = jobs.iterator();
        while(it.hasNext()) {
          Date jobDate = it.next();
          assert currentDate.after(jobDate);
          if ((currentDate.getTime() - jobDate.getTime()) > 1000 * duration) {
            // Delete jobs for the user address that have aged out.
            it.remove();
          }
        }
        if (jobs.size() == 0) {
          // All the jobs for this user address have been removed.
          // Drop this user address from the table.
          jobTable.remove(key);
        }
      }
    }
}

public class JobIPValve extends ValveBase {

    private TimeLimitedJobTable jobTable;
    private ArrayList<String> whiteList;
    private long duration = 3600; // 1 hour
    private int jobLimit = 5;

    @Override protected void initInternal() throws LifecycleException {
        jobTable = new TimeLimitedJobTable(duration, jobLimit);
        System.out.println("Initializing JobIPValve");
    }

    public long getDuration() {
      return this.duration;
    }

    public void setDuration(long duration) {
      this.duration = duration;
    }

    public int getJobLimit() {
      return this.jobLimit;
    }

    public void setJobLimit(String jobLimit) {
      this.jobLimit = Integer.parseInt(jobLimit);
    }

    public String getWhiteList() {
      final StringBuilder buff = new StringBuilder();
      for (String s : this.whiteList) {
          buff.append(s).append(" ");
      }
      return buff.toString();
    }

    public void setWhiteList(String whiteList) {
      this.whiteList = new ArrayList<String>(Arrays.asList(whiteList.split("\\s+")));
    }

    private boolean addressInWhiteList(String address) {
      for (String string : this.whiteList) {
        if(string.matches(address)){
          return true;
        }
      }
      return false;
    }

    @Override
    public void invoke(Request request, Response response) throws IOException, ServletException {
        String remoteAddress = request.getRemoteAddr();
        String requestUri = request.getRequestURI();
        String method = request.getMethod();
        if (!addressInWhiteList(remoteAddress) && method.equals("POST")) {
          if (!jobTable.addJob(remoteAddress)) {
            int numJobs = jobTable.getNumJobs(remoteAddress);
            if (numJobs >= jobLimit) {
              System.out.println("Remote host " + remoteAddress 
                + " reached the job submission rate limit."
                +  Integer.toString(jobTable.getNumJobs(remoteAddress)) 
                + " jobs in " + String.valueOf(duration) + " seconds.");
              response.sendError(
                response.SC_SERVICE_UNAVAILABLE, 
                "The server limits job submisions to "
                + String.valueOf(jobLimit) 
                + " jobs in " + String.valueOf(duration) + " seconds."
                + " Please wait at least "
                + String.valueOf(jobTable.getMinWait(remoteAddress))
                + " seconds before submitting another job."
              );
            }
            else {
              response.sendError(
                response.SC_SERVICE_UNAVAILABLE, 
                "The server has encountered a problem recording your job "
                + " with the JobIPValve component. "
                + "Please contact the server administrator with the text of this message."
              );
            }
          }
        }
        Valve nextValve = getNext();
        if(nextValve!=null){
            nextValve.invoke(request, response);
        }
   }
}
