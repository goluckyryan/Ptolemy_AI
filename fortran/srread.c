/*
 * 8/10/87 - add support for opening files append no truncation - t.h.
 * 2/24/88 - use fopen, fread, fwrite, and fclose - t.h.
 * 3/6/89 - add snwrte routine - t.h.
 * 10/6/89 - added uscores include - t.h.
 * 8/2/90 - expand srread to read 300 columns - t.h.
 * 9/6/90 - add swdlte routine - t.h.
 * 10/7/90 - add maxchr arg to srread - t.h.
 * 11/8/91 - handle lines that are not terminated by \n - t.h.
 * 1/20/93 - same as above for RS/6000 - t.h.
 * 12/3/93 - do uscores.h stuff in uscores ifdef - t.h.
 * 3/11/94 - fix up for Solaris 2.x - t.h.
 * 5/20/94 - correctly blank out newline in srread - t.h.
 *
 *  This file contains the routines used to open, close, read, and write
 * sequential files.
 * call as:
 *   srread(idcb, icode, line, numchr, maxchr, rtype)
 *   swrite(idcb, icode, line, numchr)
 *   snwrte(idcb, icode, line, numchr)
 *   sropen(idcb, icode, fname)
 *   srclse(idcb)
 *   swopen(idcb, icode, fname, oversw)
 *   swclse(idcb)
 *   swdlte(idcb, icode, fname)
 */
#ifdef system5
#include <sys/types.h>
#include <unistd.h>
#endif
#include  <sys/file.h>
#include  <errno.h>
#include  <stdio.h>
#include  <string.h>

#ifndef FALSE
#define FALSE	0
#define TRUE	1
#endif

/* ****** somehow the following is used with cc on redhat 9
****** but it should not be.  It is commented out and
****** hence AIX probably will not work
*/

/* 
#ifndef uscores
#define srread_ srread
#define swrite_ swrite
#define snwrte_ snwrte
#define sropen_ sropen
#define srclse_ srclse
#define swopen_ swopen
#define swclse_ swclse
#define swdlte_ swdlte
#endif
*/

extern int errno;

/* read numchr characters into dummy from the file pointed to by idcb */
/* in the method prescribed by the value of rtype {0 = character data;
 * 1 = binary data}  */
void srread_(idcb, icode, line, numchr, maxchr, rtype)
int *idcb, *icode, *numchr, *maxchr, *rtype;
char *line;
{
	int c, eof;
	*icode = 0;
	switch (*rtype) {
	case 0: 	/* character data */
		if (fgets(line,*maxchr,(FILE *) *idcb) == 0) *icode = 3;
		else {
			*numchr = strlen(line)==0 ? 0 : strlen(line)-1;
			if (*numchr == 0) return;
			/* if we are not at the end of a line, we have room */
			/* for one more character  */
			eof = FALSE;
			if (line[*numchr] != '\n') {	
				*numchr += 1; 		/* read an extra char */
				c = getc((FILE *) *idcb);
				if (c != EOF) {
					line[*numchr] = c;
				} else {
					eof = TRUE;
				}
			}
			if (line[*numchr] == '\n' && ! eof) {
				line[*numchr] = ' '; /* remove newline */
			} else if (! eof) {		/* move to next line */
				*numchr += 1;
				while ((c = getc((FILE *) *idcb)) != '\n' &&
					c != EOF)
					;
			}
		}
		break;
	case 1: 	/* binary data */
		if (fread(line, 80, 1, (FILE *) *idcb) == 0) {
                      *icode = 3;
                      return;
                }
                *numchr = 80;
                /* advance the pointer past the return */
                fseek( (FILE *) *idcb, 1, 1); 
		break;
	}
}

/* write numchr characters from line to the file pointed to by idcb with
 * a trailing \n */
void swrite_(idcb, icode, line, numchr)
int *idcb, *icode, *numchr;
char *line;
{
	*icode = 0;
	if (fwrite(line, 1, *numchr, (FILE *) *idcb) != *numchr) *icode = errno;
	fwrite("\n", 1, 1, (FILE *) *idcb);
}

/* write numchr characters from line to the file pointed to by idcb without
 * a trailing \n */
void snwrte_(idcb, icode, line, numchr)
int *idcb, *icode, *numchr;
char *line;
{
	*icode = 0;
	if (fwrite(line, 1, *numchr, (FILE *) *idcb) != *numchr) *icode = errno;
}

/* open the file fname for reading and put its file descriptor in idcb */
void sropen_(idcb,icode,fname)

int *idcb, *icode;
char *fname;
{
	FILE *fd;
	*icode = 0;
	if ((fd = fopen(fname,"r")) == NULL) { /* open for read only */
		*icode = errno;
		return;
	}
	*idcb = (int) fd;
}

/* close the file pointed to by idcb */
void srclse_(idcb)
int *idcb;
{
	fclose((FILE *) *idcb);
}

/* open the file fname for writing and put its file descriptor in idcb
 * make sure that all writes are appended. */
void swopen_(idcb,icode,fname,oversw,append)
int *idcb, *icode, *oversw, *append;
char *fname;
{
	FILE *fd;
	*icode = 0;

	/* if the file exists and we cannot overwrite return an error */
	if (*oversw == 0 && access(fname,F_OK) == 0) {
		*icode = EEXIST;
		return;
	}
	if (*append != 0) {   /* open for append */
		if ((fd = fopen(fname,"a")) == NULL){
			*icode = errno;   /* an error occured */
			return;
		}
		*idcb = (int) fd;
		return;
	}
	if ((fd = fopen(fname,"w")) == NULL){
		*icode = errno;    /* an error occured */
		return;
	}
	*idcb = (int) fd; /* store the file descriptor for later use */
	return;
}

/* close the file pointed to by idcb */
void swclse_(idcb)
int *idcb;
{
	fclose((FILE *) *idcb);
}


/* remove the file fname */
void swdlte_(idcb, icode, fname)
char *fname;
int *idcb, *icode;
{
	fclose((FILE *) *idcb);
	if (unlink(fname) == -1) *icode = errno;
	else *icode = 0;
	return;
}
