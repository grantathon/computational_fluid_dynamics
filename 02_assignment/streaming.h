#ifndef _STREAMING_H_
#define _STREAMING_H_

#define FLUID 0
#define NO_SLIP 1
#define MOVING_WALL 2

/** carries out the streaming step and writes the respective distribution functions from
 *  collideField to streamField.
 */
void doStreaming(double *collideField, double *streamField,int *flagField,int xlength);

#endif

