#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main()
{
	FILE * fp = NULL;
	FILE * ofp = NULL;
	char * line = NULL;
	char * subj = "Pt00";
	char * electrode = "Z00";
	char * filename = "Pt00.mni";
	char buf[9];
	int x,y,z;
	unsigned int cur_subj_num, prev_subj_num;
	size_t len = 0;
	ssize_t read;
	/* Notes to self:
	 * size_t holds an index (unsigned integer). It is used instead of uint
	 * for semantic and portability reasons.
	 * ssize_t holds an index, like size_t, but instead of being purely
	 * unsigned, it allows the value -1 as well which is used as an error
	 * code.  It's useful for storing the line length read from a file,
	 * because there may be a read error and in such cases the length will
	 * be returned as -1. */
	
	fp = fopen("MNI_basal_electrodes_Pt01_10_w_label.csv", "r");
	/* This program is designed to work only on this file, in a precise
	 * formatting. The first row is a header. There are 4 columns
	 * (subject,electrode,x,y,z). Columns are delimited with commas.
	 * Despite being MNI coordinates (which reference continuous space and
	 * so might naturally require floating point precision), the x,y,z
	 * values in this file at positive and negtive integers. Each row in
	 * the file corresponds to a different electrode. The file is already
	 * sorted by subject, so any time the value in the first column changes
	 * a new file should be opened (and the previously opened file closed,
	 * if that needs doing). The only purpose of this program is to write
	 * out the coordinates for each subject so that they can be processed
	 * with AFNI's 3dUndump to produce MRI volumes. */
	if (fp == NULL)
	{
		perror("Failed to open input file for reading");
		exit(EXIT_FAILURE);
	}

	// Skip first line (exit on error).
	if ((read = getline(&line, &len, fp)) == -1)
	{
		perror("Failed to read from input file");
		exit(EXIT_FAILURE);
	}

	// Read the file line by line.
	prev_subj_num = 0;
	while ((read = getline(&line, &len, fp)) != -1) // -1 is error code
	{
		subj = strtok(line, ",");
		electrode = strtok(NULL, ",");
		x = atoi(strtok(NULL, ","));
		y = atoi(strtok(NULL, ","));
		z = atoi(strtok(NULL, ","));
		
		/* Extract subject number from subject string. */ 
		sscanf(subj, "Pt%02d", &cur_subj_num);
			
		if (cur_subj_num > prev_subj_num)
		{
			/* If ofp is pointing to a file descriptor, close the
			 * file. */
			if (ofp != NULL)
				fclose(ofp);

			/* If sprintf cannot create the filename properly,
			 * catch and exit the program. Note that sprintf needs
			 * to assign it's output to a preallocated buffer,
			 * rather than just a pointer to a string. While many
			 * functions seem to allocate memory for you, I guess
			 * sprintf does not. */
			if ((sprintf(buf, "%s.mni", subj)) < 0)
			{
				perror("Failed to generate filename string");
				exit(EXIT_FAILURE);
			}
			/* But a pointer to the result from sprintf, stored in
			 * the buffer, can be stored in my
			 * semantically-meaningful filename variable.  */
			filename = buf;

			/* Try to open a file for writing with the filename
			 * defined in the previous step. If this failed, ofp
			 * will be NULL. Catch and exit the program. This step
			 * assumes that the rows of the input file are sorted
			 * by subject, because every time the file is opened
			 * for writing it will begin writing from the
			 * beginning, discarding anything written previously to
			 * the file. */
			ofp = fopen(filename, "w");
			if (ofp == NULL)
			{
				perror("Failed to open file for writing");
				exit(EXIT_FAILURE);
			}
			prev_subj_num = cur_subj_num;
		}
		fprintf(ofp, "%d %d %d\n", x,y,z);
	}	
	fclose(ofp);
	fclose(fp);
	/* Folks on SO mention that this bit is not strictly necessarily, but
	 * the C manual includes the free(line) bit so I just extrapolated that to all
	 * my strings that are stored on the heap. Note that the buf is
	 * allocated locally and so is not (and cannot) be freed. It is
	 * naturally cleaned up when the program exits.*/
	// STRIKE THAT. Freeing other string lest to a crash. I am going to go
	// ahead and not free anything.
	// if (line)
	//	free(line);
	exit(EXIT_SUCCESS);
}
