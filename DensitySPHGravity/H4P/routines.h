/* 
** routines.h
**
** Header file for routines.c
*/

void initialise_parameters(SI (*));
void initialise_black_hole(PARTICLE (*));
void check_main_parameters(const SI (*));
void calculate_parameters(const GI (*), SI (*)); 
void initialise_gridr(GI (*), PARTICLE (*), SI (*));
void initialise_griddf(const GI (*), SI (*));
void initialise_structure(const GI (*), SI (*));
void set_positions(SI (*)); 
void set_velocities(const GI (*), SI (*));
void set_attributes(const GI (*), SI (*));
void double_particles(SI (*));
void calculate_stuff(GI (*), PARTICLE (*), SI (*));
void write_griddf(SI (*), FILE (*));
void write_gridr(GRIDR (*), FILE (*));
void usage(void);
