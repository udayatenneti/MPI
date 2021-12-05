
#ifndef _BODY_H_N_BODY_
#define _BODY_H_N_BODY_

typedef struct {
    double pos[2];
    double vel[2];
    double m;
    double work;
    int    idx;
} Body;

#endif // _BODY_H_N_BODY_
