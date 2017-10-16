/* ==================================================================================
  Program to count the number of independent members in large ensemble

  Language: C++

  Modification history
  --------------------
  09/09/17 New Code. Ross Bannister, r.n.bannister@reading.ac.uk

  On Ubuntu machines
  g++ -I/usr/include IndependentMemberCount.cpp -L/usr/lib -lnetcdf_c++ -lnetcdf

  On Reading university linux
  g++ -I/opt/graphics/64/include IndependentMemberCount.cpp -L/opt/graphics/64/lib -lnetcdf_c++ -lnetcdf-4.0 -lhdf5_hl -lhdf5

  To run
  ./a.out

  GNU GENERAL PUBLIC LICENSE
  Version 3, 29 June 2007
  See file "LICENSE" for details

  =============================================================================== */

  #include <stdio.h>
  #include <math.h>
  #include <stdlib.h>
  #include <string.h>
  #include <netcdf.h>

/* -------------------------------------------------------------------------------
   Global constants
   ------------------------------------------------------------------------------- */
  const int    lenfilen       = 256;
  const int    Nx             = 360;
  const int    Ny             = 288;
  const int    nens           = 93;     // Number of members in ensemble
  const int    ens_type       = 1;      // 1 = large ensemble
                                        // 2 = ETKF ensemble
  const bool   pert_from_mean = false;  // Take perts from mean, or from control
  const bool   check_results  = false;

/* -------------------------------------------------------------------------------
   Global variables
   ------------------------------------------------------------------------------- */

  char   param      = 't';  // t or q
  int    lev        = 35;   // Level of interest for temperature
//  char   param      = 'q';  // t or q
//  int    lev        = 10;   // Level of interest for specific humidity


/* -------------------------------------------------------------------------------
   Define the important data types
   ------------------------------------------------------------------------------- */
  struct field_type
  { double xvals[::Nx];
    double yvals[::Ny];
    double value[::Nx][::Ny];
  };

  struct GramSchmidt_type
  { double alpha[::nens][::nens];
    double N[::nens];
    double N_hat[::nens];
  };


/* -------------------------------------------------------------------------------
   Function declarations
   ------------------------------------------------------------------------------- */

double InnerProd ( struct field_type *x1,
                   struct field_type *x2 );

void Divide ( struct field_type *x,
              double            N );

void Copy ( struct field_type *x1,
            struct field_type *x2 );

void Intit_alpha ( struct GramSchmidt_type *GS);

void Combine ( struct field_type *x1,
               struct field_type *x2,
               double            fac );

void ComputePerts ( struct field_type x[],
        	    int               nmems,
        	    bool              pert_from_mean,
        	    int               *nperts );

void ReadEns (
       int    ens_type,                // 1 = large ens, 2=etkf ens
       char   path[],                  // Pathname
       char   param,                   // t or q
       int    time_counter,            // time index (hour)
       int    lev,                     // level index
       struct field_type ens[::nens] );// out ensemble data



/* -------------------------------------------------------------------------------
   Main part of the program
   ------------------------------------------------------------------------------- */
int main ()

{ // Declare variables
  char                     inpath[::lenfilen];
  char                     outfile[::lenfilen];
  char                     pert_or_cntl[8];
  struct field_type        *ens = NULL;
  struct field_type        *ortho = NULL;
  struct field_type        intermediate;
  struct GramSchmidt_type  GS;
  double                   value;
  int                      nperts;
  int                      i, j, th, tm;
  FILE*                    op = NULL;
  int                      t_hr_min, t_hr_max, t_mn_min, t_mn_max, subhr_intervals;
  bool                     go;
  int                      time_counter = -1;

  // Allocate memory from the heap for the ensemble and the set of orthogonal vectors
  ens   = new struct field_type[::nens];
  ortho = new struct field_type[::nens];

  // Determine part of filename to indicate pert from mean or from control
  if (::pert_from_mean)
  { sprintf (pert_or_cntl, "frommean");
  }
  else
  { sprintf (pert_or_cntl, "fromcntl");
  }

  // Set-up some stuff depending on ensemble type
  switch (::ens_type)
  { case 1: // large ensemble
      t_hr_min        = 8;
      t_mn_min        = 0;
      t_hr_max        = 17;
      t_mn_max        = 0;
      subhr_intervals = 0;
      sprintf (inpath, "/export/diamet/raid3/stefano/cases/largens");
      break;
    case 2: // ETKF ensemble
      t_hr_min        = 11;
      t_mn_min        = 30;
      t_hr_max        = 15;
      t_mn_max        = 30;
      subhr_intervals = 1;
      break;
  }

  // Loop over time
  for (th=t_hr_min; th<=t_hr_max; th++)
  { for (tm=0; tm<=subhr_intervals; tm++)
    { // Do we deal with this ensemble?
      go = true;
      if ((th==t_hr_min) && (tm*30<t_mn_min))
      { go = false;
      }
      if ((th==t_hr_max) && (tm*30>t_mn_max))
      { go = false;
      }
      // printf ("th = %02d, tm = %02d, go = %c\n", th, tm, go);
      if (go)
      { // Set alphas to zero
        Intit_alpha ( &GS );

        time_counter++;
        // Read-in the full ensemble
        printf ("Reading in the ensemble\n");
        if (ens_type == 2)
        { sprintf (inpath, "/export/diamet/raid3/ross/Ensembles/2011-09-20-%02u%02u_24mems_NoRP",
                           th, tm*30);
        }

        // Read-in the ensmble
        ReadEns ( ::ens_type,
                  inpath,
                  ::param,
                  time_counter,
                  ::lev,
                  ens );

        // Compute the perturbations
        ComputePerts ( ens,
        	       ::nens,
        	       ::pert_from_mean,
        	       &nperts );
        printf ("There are %u perturbations\n", nperts);

        // Determine the output file name, and open
        switch (::ens_type)
        { case 1:
            sprintf (outfile, "/export/diamet/raid3/ross/tmp/Normalizations_large_%s_%02u%02u_%c",
            pert_or_cntl, th, tm*30, ::param);
            break;
          case 2:
            sprintf (outfile, "/export/diamet/raid3/ross/tmp/Normalizations_etkf_%s_%02u%02u_%c",
            pert_or_cntl, th, tm*30, ::param);
            break;
        }
        op = fopen (outfile, "w");

        // Set-up the first vector (vector 0)
        printf ("Setting up orthogonal vectors\n");
        GS.N[0] = sqrt(InnerProd ( &(ens[0]),
                                   &(ens[0]) ));
        GS.N_hat[0] = 1.0;

        Copy   ( &(ens[0]),
                 &(ortho[0]) );
        Divide ( &(ortho[0]),
                 GS.N[0] );

        fprintf (op, "%u  %f\n", 0, 1.0);

        // Find the later vectors
        for (i=1; i<nperts; i++)
        { // Find the ith vector

          // Normalize the ith member, and then place into the ortho structure for now
          GS.N[i] = sqrt(InnerProd ( &(ens[i]),
                                     &(ens[i]) ));
          Copy   ( &(ens[i]),
                   &(ortho[i]) );
          Divide ( &(ortho[i]),
                   GS.N[i] );

          // Find the alphas for this vector
          for (j=0; j<i; j++)
          { GS.alpha[i][j] = InnerProd ( &(ortho[i]),  // actually ens[i]/GS.N[i]
                                         &(ortho[j]) );
          }

          // Build the rest of the ith orthogonal vector
          for (j=0; j<i; j++)
          { Combine ( &(ortho[i]),
                      &(ortho[j]),
                      GS.alpha[i][j] );
          }

          // x_hat[i] contains the un-normalized orthoginal vector, so normalize
          GS.N_hat[i] = sqrt(InnerProd ( &(ortho[i]),
                                         &(ortho[i]) ));
          Divide ( &(ortho[i]),
                   GS.N_hat[i] );

          fprintf (op, "%u  %f\n", i, GS.N_hat[i]);

        }

        fclose (op);

        if (check_results)
        { // Check to see if ensemble is orthonormal
          for (i=0; i<nperts; i++)
          { printf ("%02u :", i);
            for (j=0; j<nperts; j++)
            { printf ("  %f", InnerProd ( &(ortho[i]),
                                          &(ortho[j]) ) );
            }
            printf ("\n");
          }
        }
      }
    }
  }

  // Deallocate the data structures
  delete[] ens, ortho;

}






/* -------------------------------------------------------------------------------
   Compute the inner product between two fields
   ------------------------------------------------------------------------------- */
double InnerProd ( struct field_type *x1,
                   struct field_type *x2 )
{ // Declare variables
  int    i, j;
  double sum;

  sum = 0.0;
  for (i=0; i<(::Nx); i++)
  { for (j=0; j<(::Ny); j++)
    { sum += (*x1).value[i][j] * (*x2).value[i][j];
    }
  }

  return sum;
}



/* -------------------------------------------------------------------------------
   Divide a field by a number
   ------------------------------------------------------------------------------- */
void Divide ( struct field_type *x,
              double            N )
{ // Declare variables
  int    i, j;
  double fac;

  fac = 1.0 / N;
  for (i=0; i<(::Nx); i++)
  { for (j=0; j<(::Ny); j++)
    { (*x).value[i][j] *= fac;
    }
  }
}



/* -------------------------------------------------------------------------------
   Copy a field
   ------------------------------------------------------------------------------- */
void Copy ( struct field_type *x1,
            struct field_type *x2 )
{ // Declare variables
  int    i, j;

  for (i=0; i<(::Nx); i++)
  { for (j=0; j<(::Ny); j++)
    { (*x2).value[i][j] = (*x1).value[i][j];
    }
  }

  for (i=0; i<(::Nx); i++)
  { (*x2).xvals[i] = (*x1).xvals[i];
  }
  for (j=0; j<(::Ny); j++)
  { (*x2).yvals[j] = (*x1).yvals[j];
  }

}


/* -------------------------------------------------------------------------------
   Set alphas to zero
   ------------------------------------------------------------------------------- */
void Intit_alpha ( struct GramSchmidt_type *GS)
{ // Declare variables
  int    i, j;

  for (i=0; i<(::nens); i++)
  { for (j=0; j<(::nens); j++)
    { (*GS).alpha[i][j] = 0.0;
    }
  }

}


/* -------------------------------------------------------------------------------
   Combine fields as follows x1 := x1 - fac * x2
   ------------------------------------------------------------------------------- */
void Combine ( struct field_type *x1,
               struct field_type *x2,
               double            fac )
{ // Declare variables
  int    i, j;

  for (i=0; i<(::Nx); i++)
  { for (j=0; j<(::Ny); j++)
    { (*x1).value[i][j] -= fac * (*x2).value[i][j];
    }
  }

}


/* -------------------------------------------------------------------------------
   Compute perturbations
   ------------------------------------------------------------------------------- */
void ComputePerts ( struct field_type x[],
        	    int               nmems,
        	    bool              pert_from_mean,
        	    int               *nperts )
{ // Declare variables
  int               i, j, member, pert, offset, pertoffset;
  double            mean, fac;
  struct field_type ref;

  if (pert_from_mean)
  { // Compute perturbations from the mean, so compute the mean
    printf ("Computing perturbations from the mean\n");
    *nperts = nmems;
    offset    = 0;
    fac = 1.0 / double(nmems);
    for (i=0; i<(::Nx); i++)
    { for (j=0; j<(::Ny); j++)
      { mean = 0.0;
        for (member=0; member<nmems; member++)
        { mean += x[member].value[i][j];
        }
        ref.value[i][j] = fac * mean;
      }
    }
  }
  else
  { // The reference state is the control
    printf ("Computing perturbations from the control member\n");
    *nperts = nmems - 1;
    offset    = 1;
    Copy ( &(x[0]),
           &ref );
  }

  // Compute the perts from the reference state
  for (pert=0; pert<(*nperts); pert++)
  { pertoffset = pert + offset;
    for (i=0; i<(::Nx); i++)
    { for (j=0; j<(::Ny); j++)
      { x[pert].value[i][j] = x[pertoffset].value[i][j] - ref.value[i][j];
      }
    }
  }


}


/* -------------------------------------------------------------------------------
   Read ensemble members
   ------------------------------------------------------------------------------- */
void ReadEns (
       int    ens_type,                // 1 = large ens, 2=etkf ens
       char   path[],                  // Pathname
       char   param,                   // t or q
       int    time_counter,            // time counter
       int    lev,                     // level index
       struct field_type ens[::nens] ) // out ensemble data

{ // Declare local variables
  char              filename[::lenfilen];
  int               ensloop, x, y;
  int               ncid, ierr, ierr1, ierr2, ierr3;
  int               varidtheta, varidexner, varidq;
  struct field_type *theta = NULL;
  struct field_type *exner = NULL;
  size_t            start[4], count[4];
  float             lineofdata[::Ny];


  if (param == 't')
  { theta = new struct field_type;
    exner = new struct field_type;
  }

  // Set bounds for retrieving the variables
                    count[3] = 1;         // x
  start[2] = 0;     count[2] = ::Ny;      // y
  start[1] = lev;   count[1] = 1;         // z
                    count[0] = 1;         // t

  switch (ens_type)
  { case 1:
      start[0] = time_counter;
      break;
    case 2:
      start[0] = 0;
      break;
  }


  for (ensloop=0; ensloop<(::nens); ensloop++)
  { // Construct the filename
    switch (ens_type)
    { case 1:
        sprintf (filename, "%s/qwq107.oper%02d.pp1.full.nc", path, ensloop);
        break;
      case 2:
        sprintf (filename, "%s/Member%02d.noRP.nc", path, ensloop);
        break;
    }
    printf ("Time index %u, Reading file %s\n", time_counter, filename);
    ierr = nc_open (filename,
                    NC_NOWRITE,
                    &ncid);
    if (ierr != 0)
    { printf ("Error opening file %s\n%s\nExiting\n", filename, nc_strerror(ierr));
      exit(0);
    }

    // Get variable ids
    if (param == 't')
    { ierr1 = nc_inq_varid (ncid,
                            "theta",
                            &varidtheta);
      ierr2 = nc_inq_varid (ncid,
                            "field7",
                            &varidexner);

      // Read-in temperature
      for (x=0; x<(::Nx); x++)
      { start[3] = x;
        ierr     = nc_get_vara_float (ncid,
                                      varidtheta,
                                      start,
                                      count,
                                      lineofdata);
        if (ierr != 0)
        { printf ("Error reading theta\n%s\nExiting\n", nc_strerror(ierr));
          exit(0);
        }
        // Insert these data into Field array
        for (y=0; y<(::Ny); y++)
        { (*theta).value[x][y] = lineofdata[y];
        }
      }

      // Read-in exner
      for (x=0; x<(::Nx); x++)
      { start[3] = x;
        ierr     = nc_get_vara_float (ncid,
                                      varidexner,
                                      start,
                                      count,
                                      lineofdata);
        if (ierr != 0)
        { printf ("Error reading exner\n%s\nExiting\n", nc_strerror(ierr));
          exit(0);
        }
        // Insert these data into Field array
        for (y=0; y<(::Ny); y++)
        { (*exner).value[x][y] = lineofdata[y];
        }
      }

      // Compute temperature
      for (x=0; x<(::Nx); x++)
      { for (y=0; y<(::Ny); y++)
        { ens[ensloop].value[x][y] = (*theta).value[x][y] * (*exner).value[x][y];
        }
      }
    }


    if (param == 'q')
    { ierr1 = nc_inq_varid (ncid,
                            "q",
                            &varidq);
      // Read-in specific humidity
      for (x=0; x<(::Nx); x++)
      { start[3] = x;
        ierr     = nc_get_vara_float (ncid,
                                      varidq,
                                      start,
                                      count,
                                      lineofdata);
        if (ierr != 0)
        { printf ("Error reading q\n%s\nExiting\n", nc_strerror(ierr));
          exit(0);
        }
        // Insert these data into Field array
        for (y=0; y<(::Ny); y++)
        { ens[ensloop].value[x][y] = 1000.0 * lineofdata[y];
        }
      }
    }
  }

  // Tidy up
  if (param == 't')
  { delete theta;
    delete exner;
  }
}
