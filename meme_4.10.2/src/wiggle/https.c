/* Connect via https. */

#ifdef USE_SSL

#include "openssl/ssl.h"
#include "openssl/err.h"

#include <sys/socket.h>
#include <unistd.h>
#include <signal.h>

#include "common.h"
#include "errabort.h"


static void xerrno(char *msg)
{
fprintf(stderr, "%s : %s\n", strerror(errno), msg); fflush(stderr);
}

static void xerr(char *msg)
{
fprintf(stderr, "%s\n", msg); fflush(stderr);
}


int netConnectHttps(char *hostName, int port)
/* Return socket for https connection with server or -1 if error. */
{

fflush(stdin);
fflush(stdout);
fflush(stderr);

int sv[2]; /* the pair of socket descriptors */

socketpair(AF_UNIX, SOCK_STREAM, 0, sv);

int pid = fork();

if (pid < 0)
    errnoAbort("can't fork in netConnectHttps");
if (pid == 0)
    {
    /* child */

    fclose(stdin);
    fclose(stdout);

    close(sv[0]);  /* close unused half of pipe */

    /* close other file descriptors */
    int fd=0;
    for (fd = STDERR_FILENO+1; fd < 64; fd++)
      if (fd != sv[1])
  	close(fd);

    char hostnameProto[256];

    BIO *sbio;
    SSL_CTX *ctx;
    SSL *ssl;

    SSL_library_init();

    ERR_load_crypto_strings();
    ERR_load_SSL_strings();
    OpenSSL_add_all_algorithms();

    /* We would seed the PRNG here if the platform didn't
    * do it automatically
    */

    ctx = SSL_CTX_new(SSLv23_client_method());

    fd_set readfds;
    fd_set writefds;
    int err;
    struct timeval tv;


    /* future extension: checking certificates 

    char *certFile = NULL;
    char *certPath = NULL;
    if (certFile || certPath)
	{
	SSL_CTX_load_verify_locations(ctx,certFile,certPath);
    #if (OPENSSL_VERSION_NUMBER < 0x0090600fL)
	SSL_CTX_set_verify_depth(ctx,1);
    #endif
	}

    */

    /* We'd normally set some stuff like the verify paths and
    * mode here because as things stand this will connect to
    * any server whose certificate is signed by any CA.
    */

    sbio = BIO_new_ssl_connect(ctx);

    BIO_get_ssl(sbio, &ssl);
    if(!ssl) 
	{
	xerr("Can't locate SSL pointer");
	goto cleanup;
	}

    /* Don't want any retries since we are non-blocking bio now */
    //SSL_set_mode(ssl, SSL_MODE_AUTO_RETRY);

    /* We might want to do other things with ssl here */

    safef(hostnameProto,sizeof(hostnameProto),"%s:%d",hostName,port);
    BIO_set_conn_hostname(sbio, hostnameProto);

    BIO_set_nbio(sbio, 1);     /* non-blocking mode */


    while (1) 
	{
	if (BIO_do_connect(sbio) == 1) 
	    {
	    break;  /* Connected */
	    }
	if (! BIO_should_retry(sbio)) 
	    {
	    xerr("BIO_do_connect() failed");
	    char s[256];	
	    safef(s, sizeof s, "SSL error: %s", ERR_reason_error_string(ERR_get_error()));
	    xerr(s);
	    goto cleanup;
	    }

	fd = BIO_get_fd(sbio, NULL);
	if (fd == -1) 
	    {
	    xerr("unable to get BIO descriptor");
	    goto cleanup;
	    }

	FD_ZERO(&readfds);
	FD_ZERO(&writefds);
	if (BIO_should_read(sbio)) 
	    {
	    FD_SET(fd, &readfds);
	    }
	else if (BIO_should_write(sbio)) 
	    {
	    FD_SET(fd, &writefds);
	    }
	else 
	    {  /* BIO_should_io_special() */
	    FD_SET(fd, &readfds);
	    FD_SET(fd, &writefds);
	    }
	tv.tv_sec = 10;  // timeout
	tv.tv_usec = 0;

	err = select(fd + 1, &readfds, &writefds, NULL, &tv);

	if (err < 0) 
	    {
	    xerr("select() error");
	    goto cleanup;
	    }

	if (err == 0) 
	    {
	    char s[256];	
	    safef(s, sizeof s, "connection timeout to %s", hostName);
	    xerr(s);
	    goto cleanup;
	    }
	}



    /* future extension: checking certificates 

    if (certFile || certPath)
	if (!check_cert(ssl, host))
	    return -1;

    */

    /* Could examine ssl here to get connection info */

    /* we need to wait on both the user's socket and the BIO SSL socket 
     * to see if we need to ferry data from one to the other */


    char sbuf[32768];  // socket buffer sv[1] to user
    char bbuf[32768];  // bio buffer
    int srd = 0;
    int swt = 0;
    int brd = 0;
    int bwt = 0;
    while (1) 
	{

	// Do NOT move this outside the while loop. 
	/* Get underlying file descriptor, needed for select call */
	fd = BIO_get_fd(sbio, NULL);
	if (fd == -1) 
	    {
	    xerr("BIO doesn't seem to be initialized in https, unable to get descriptor.");
	    goto cleanup;
	    }


	FD_ZERO(&readfds);
	FD_ZERO(&writefds);

	if (brd == 0)
    	    FD_SET(fd, &readfds);
	if (swt < srd)
	    FD_SET(fd, &writefds);
	if (srd == 0)
    	    FD_SET(sv[1], &readfds);

	tv.tv_sec = 10;   // timeout
	tv.tv_usec = 0;

	err = select(max(fd,sv[1]) + 1, &readfds, &writefds, NULL, &tv);
    
	/* Evaluate select() return code */
	if (err < 0) 
	    {
	    xerr("error during select()");
	    goto cleanup;
	    }
	else if (err == 0) 
	    {
	    /* Timed out - just quit */
	    xerr("https timeout expired");
	    goto cleanup;
	    }

	else 
	    {
	    if (FD_ISSET(sv[1], &readfds))
		{
		swt = 0;
		srd = read(sv[1], sbuf, 32768);
		if (srd == -1)
		    {
		    if (errno != 104) // udcCache often closes causing "Connection reset by peer"
			xerrno("error reading https socket");
		    goto cleanup;
		    }
		if (srd == 0) 
		    break;  // user closed socket, we are done
		}

	    if (FD_ISSET(fd, &writefds))
		{
		int swtx = BIO_write(sbio, sbuf+swt, srd-swt);
		if (swtx <= 0)
		    {
		    if (!BIO_should_write(sbio))
			{
			ERR_print_errors_fp(stderr);
			xerr("Error writing SSL connection");
			goto cleanup;
			}
		    }
		else
		    {
		    swt += swtx;
		    if (swt >= srd)
			{
			swt = 0;
			srd = 0;
			}
		    }
		}

	    if (FD_ISSET(fd, &readfds))
		{
		bwt = 0;
		brd = BIO_read(sbio, bbuf, 32768);

		if (brd <= 0)
		    {
		    if (BIO_should_read(sbio))
			{
			brd = 0;
			continue;
			}
		    else
			{
			if (brd == 0) break;
			ERR_print_errors_fp(stderr);
			xerr("Error reading SSL connection");
			goto cleanup;
			}
		    }
                // write the https data received immediately back on socket to user, and it's ok if it blocks.
		while(bwt < brd)
		    {
		    int bwtx = write(sv[1], bbuf+bwt, brd-bwt);
		    if (bwtx == -1)
			{
			if ((errno != 104)  // udcCache often closes causing "Connection reset by peer"
			 && (errno !=  32)) // udcCache often closes causing "Broken pipe"
			    xerrno("error writing https data back to user socket");
			goto cleanup;
			}
		    bwt += bwtx;
		    }
		brd = 0;
		bwt = 0;
		}
	    }
	}

cleanup:
    BIO_free_all(sbio);
    close(sv[1]);  /* we are done with it */

    exit(0);

    /* child will never get to here */
    }

/* parent */

close(sv[1]);  /* close unused half of socket */

return sv[0];

}

#else

#include <stdarg.h>
#include "common.h"
#include "errabort.h"

int netConnectHttps(char *hostName, int port)
/* Start https connection with server or die. */
{
errAbort("No openssl available in netConnectHttps for %s : %d", hostName, port);
return -1;   /* will never get to here, make compiler happy */
}

#endif
