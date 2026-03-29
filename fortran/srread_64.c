/*
 * 64-bit compatible version of srread.c for Ptolemy
 *
 * The original srread.c stores FILE* pointers in int variables,
 * which truncates the pointer on 64-bit systems. This version
 * uses intptr_t to hold the full pointer value.
 *
 * Original by T.H. (1987-1994), 64-bit fix by Dudu (2026).
 */
#include <unistd.h>
#include <sys/file.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#ifndef FALSE
#define FALSE 0
#define TRUE  1
#endif

extern int errno;

void srread_(intptr_t *idcb, int *icode, char *line, int *numchr,
             int *maxchr, int *rtype)
{
    int c, eof;
    *icode = 0;
    switch (*rtype) {
    case 0:
        if (fgets(line, *maxchr, (FILE *) *idcb) == 0) *icode = 3;
        else {
            *numchr = strlen(line)==0 ? 0 : strlen(line)-1;
            if (*numchr == 0) return;
            eof = FALSE;
            if (line[*numchr] != '\n') {
                *numchr += 1;
                c = getc((FILE *) *idcb);
                if (c != EOF) { line[*numchr] = c; } else { eof = TRUE; }
            }
            if (line[*numchr] == '\n' && !eof) {
                line[*numchr] = ' ';
            } else if (!eof) {
                *numchr += 1;
                while ((c = getc((FILE *) *idcb)) != '\n' && c != EOF) ;
            }
        }
        break;
    case 1:
        if (fread(line, 80, 1, (FILE *) *idcb) == 0) {
            *icode = 3; return;
        }
        *numchr = 80;
        fseek((FILE *) *idcb, 1, 1);
        break;
    }
}

void swrite_(intptr_t *idcb, int *icode, char *line, int *numchr)
{
    *icode = 0;
    if ((int)fwrite(line, 1, *numchr, (FILE *) *idcb) != *numchr)
        *icode = errno;
    fwrite("\n", 1, 1, (FILE *) *idcb);
}

void snwrte_(intptr_t *idcb, int *icode, char *line, int *numchr)
{
    *icode = 0;
    if ((int)fwrite(line, 1, *numchr, (FILE *) *idcb) != *numchr)
        *icode = errno;
}

void sropen_(intptr_t *idcb, int *icode, char *fname)
{
    FILE *fd;
    *icode = 0;
    if ((fd = fopen(fname, "r")) == NULL) { *icode = errno; return; }
    *idcb = (intptr_t) fd;
}

void srclse_(intptr_t *idcb)
{
    fclose((FILE *) *idcb);
}

void swopen_(intptr_t *idcb, int *icode, char *fname, int *oversw,
             int *append)
{
    FILE *fd;
    *icode = 0;
    if (*oversw == 0 && access(fname, F_OK) == 0) {
        *icode = EEXIST; return;
    }
    if (*append != 0) {
        if ((fd = fopen(fname, "a")) == NULL) { *icode = errno; return; }
        *idcb = (intptr_t) fd;
        return;
    }
    if ((fd = fopen(fname, "w")) == NULL) { *icode = errno; return; }
    *idcb = (intptr_t) fd;
}

void swclse_(intptr_t *idcb)
{
    fclose((FILE *) *idcb);
}

void swdlte_(intptr_t *idcb, int *icode, char *fname)
{
    fclose((FILE *) *idcb);
    if (unlink(fname) == -1) *icode = errno;
    else *icode = 0;
}
