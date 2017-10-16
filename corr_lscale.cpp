/* ==================================================================================
  Program to compute lengthscales in 2D correlation data

  Language: C++

  Modification history
  --------------------
  18/11/16 New Code. Ross Bannister, r.n.bannister@reading.ac.uk

  On Ubuntu machines
  g++ -I/usr/include corr_lscale.cpp -L/usr/lib -lnetcdf_c++ -lnetcdf

  On Reading university linux
  g++ -I/opt/graphics/64/include corr_lscale.cpp -L/opt/graphics/64/lib -lnetcdf_c++ -lnetcdf-4.0 -lhdf5_hl -lhdf5

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
  const int    lenfilen   = 256;
  const int    Nx         = 360;
  const int    Ny         = 288;
  const int    nens       = 93;   // Number of members in large ensemble
  const double deltax     = 1.5;  // Grid size in km
  const int    maxiters   = 10;   // Maximum number of iterations
  const double critdiff   = 0.01; // Stopping criterion for iterations
  const int    radius     = 30;   // correlations computed at this radius
  const int    Nparams    = 7;    // Number of parameters describing the correlation function
  const double y_stddev   = 1.0;  // Defines tightness of fit to data


/* -------------------------------------------------------------------------------
   Global variables
   ------------------------------------------------------------------------------- */

  char   inpath[]   = "/export/diamet/raid3/stefano/cases/largens";
  char   outpath[]  = "/export/diamet/raid3/ross/cases/largens/raw_ross";
  int    nens_red   = 24;   // Number of members in reduced ensemble
  int    t          = 7;    // Time of interest
  //int    lev        = 35;   // Level of interest
  //char   param      = 't';  // t or q
  int    lev        = 10;   // Level of interest
  char   param      = 'q';  // t or q
  bool   choosemembers24[::nens] = {true,false,false,false,false, true,false, true,false,false,false,false,
                                    true, true,false,false,false,false,false,false,false,false, true,false,
                                    true,false,false,false,false,false, true,false,false, true,false, true,
                                   false, true,false,false,false,false, true,false, true,false,false,false,
                                   false,false,false,false, true,false,false,false,false, true,false,false,
                                    true,false, true,false,false,false,false,false,false,false, true, true,
                                    true,false,false,false,false,false, true,false,false,false, true,false,
                                   false,false,false,false,false,false, true, true,false};
  bool   choosemembers47[::nens] = { true, true, true,false, true,false, true, true, true, true, true,false,
                                    false,false,false,false, true, true,false,false,false,false,false,false,
                                     true,false, true, true, true,false, true,false, true,false, true, true,
                                    false,false, true,false, true,false,false,false,false, true, true,false,
                                     true, true, true,false, true,false, true, true, true,false, true,false,
                                    false,false, true,false,false,false,false, true, true, true, true,false,
                                    false, true,false, true, true,false, true,false,false,false,false,false,
                                     true, true,false, true, true, true, true,false, true};
  bool   choosemembers93[::nens] = { true, true, true, true, true, true, true, true, true, true, true, true,
                                     true, true, true, true, true, true, true, true, true, true, true, true,
                                     true, true, true, true, true, true, true, true, true, true, true, true,
                                     true, true, true, true, true, true, true, true, true, true, true, true,
                                     true, true, true, true, true, true, true, true, true, true, true, true,
                                     true, true, true, true, true, true, true, true, true, true, true, true,
                                     true, true, true, true, true, true, true, true, true, true, true, true,
                                     true, true, true, true, true, true, true, true, true};



/* -------------------------------------------------------------------------------
   Define the important data types
   ------------------------------------------------------------------------------- */
  struct field_type
  { double xvals[::Nx];
    double yvals[::Ny];
    double value[::Nx][::Ny];
  };

  struct sub_field_type
  { double xvals[2*::radius+1];
    double yvals[2*::radius+1];
    double value[2*::radius+1][2*::radius+1];
  };

/* -------------------------------------------------------------------------------
   Function declarations
   ------------------------------------------------------------------------------- */

void ReadEns (
       char   path[],                  // Pathname
       char   param,                   // t or q
       int    t,                       // time index
       int    lev,                     // level index
       struct field_type ens[::nens] );// out ensemble data

void writefield ( char              path[],
                  char              file[],
                  struct field_type *field );

void writesubfield ( char                  path[],
                     char                  file[],
                     struct sub_field_type *field );

void mean_stddev_calc ( struct field_type ens[],
                        int               nens,
                        struct field_type *mean,
                        struct field_type *std,
                        bool              choosemembers[] );

void corr_calc ( struct field_type     ens[],
                 int                   nens,
                 struct field_type     *std,
                 struct field_type     *mean,
                 struct sub_field_type *corr,
                 int                   x0,     // The source point
                 int                   y0,
                 int                   x1,     // Lower point
                 int                   y1,
                 int                   x2,     // Upper point
                 int                   y2,
                 bool                  choosemembers[] );

void init_field ( struct sub_field_type *field );

double exponential( double x,
                    double y,
                    double p[7] );

void fn ( struct sub_field_type *field,
          double                p[7] );

void Fmatrix ( double xvals[2*::radius+1],
               double yvals[2*::radius+1],
               double p0[7],
               double delta[7],
               double F[(2*::radius+1)*(2*::radius+1)][7] );

void FTF ( double F[2*::radius+1][7],
           double recip_y_stddev2,
           double product[7][7] );

void GaussEl ( double A[7][7],       // Matrix (modified on output)
               double x[7],          // Solution (output)
               double y[7] );        // RHS (modified on output)


void LeastSq ( double p[7],                // Parameters to optimise
               double apriori[7],          // a-priori parameter values
               double delta[7],            // Small perts to parameters (for finite differencing)
               double p_stddev[7],         // sqroot of diagonal of P
               double y_stddev,            // sqroot of diagonal of R
               double y[(2*::radius+1)*(2*::radius+1)],     // Data to fit to
               double xvals[2*::radius+1], // x-values (co-ordinates)
               double yvals[2*::radius+1], // y-values (co-ordinates)
               int    maxiters,            // Max number of iterations
               double critdiff );          // Stopping criteria




/* -------------------------------------------------------------------------------
   Main part of the program
   ------------------------------------------------------------------------------- */
int main ()

{ // Declare variables
  char                  outfile[::lenfilen];
  struct field_type     *ens = NULL;
  struct field_type     mean_full;
  struct field_type     std_full;
  struct field_type     mean_red;
  struct field_type     std_red;
  struct field_type     lengths;
  struct sub_field_type corr;
  double                corr1[(2*::radius+1)*(2*::radius+1)];
  int                   x0, y0, x_low, y_low, x_hi, y_hi, c, x, y, pp;
  int                   mem;
  double                p[7], delta[7], apriori[7], p_stddev[7];
  bool                  *choosemembers;


  // Define small parameter increments (for finite differences)
  delta[0]    = 0.001;   // x0  (centre of ellipse)
  delta[1]    = 0.001;   // y0  (centre of ellipse)
  delta[2]    = 0.001;   // Lxp (lengthscale along principle axis 1)
  delta[3]    = 0.001;   // Lyp (lengthscale along principle axis 2)
  delta[4]    = 0.001;   // C0  (amplitude)
  delta[5]    = 0.001;   // D0  (shift)
  delta[6]    = 0.001;   // rot (rotation)

  // Define the a-priori p-value
  apriori[0]  = 0.0;     // x0  (centre of ellipse)
  apriori[1]  = 0.0;     // y0  (centre of ellipse)
  apriori[2]  = 50.0;    // Lxp (lengthscale along principle axis 1)
  apriori[3]  = 50.0;    // Lyp (lengthscale along principle axis 2)
  apriori[4]  = 1.0;     // C0  (amplitude)
  apriori[5]  = 0.0;     // D0  (shift)
  apriori[6]  = 0.0;     // rot (rotation)

  // Define the a-priori standard deviation
  p_stddev[0] = 1.0;     // x0  (centre of ellipse)
  p_stddev[1] = 1.0;     // y0  (centre of ellipse)
  p_stddev[2] = 30.0;    // Lxp (lengthscale along principle axis 1)
  p_stddev[3] = 30.0;    // Lyp (lengthscale along principle axis 2)
  p_stddev[4] = 0.1;     // C0  (amplitude)
  p_stddev[5] = 0.1;     // D0  (shift)
  p_stddev[6] = 4.0;     // rot (rotation)

  // Allocate memory from the heap for the ensemble
  ens = new struct field_type[::nens];


  // Read-in the ensemble
  printf ("Reading in the full ensemble\n");
  ReadEns ( ::inpath,
            ::param,
            ::t,
            ::lev,
            ens );

  // Choose the members that are to be used for the reduced ensemble
  switch (::nens_red)
  { case 24:
      choosemembers = ::choosemembers24;
      break;
    case 47:
      choosemembers = ::choosemembers47;
      break;
    case 93:
      choosemembers = ::choosemembers93;
      break;
  }

  // Initialise the x and y values associated with the ensemble
  for (mem=0; mem<(::nens); mem++)
  { for (x0=0; x0<(::Nx); x0++)
    { ens[mem].xvals[x0] = ::deltax * float(x0);
    }
    for (y0=0; y0<(::Ny); y0++)
    { ens[mem].yvals[y0] = ::deltax * float(y0);
    }
  }


  // Output a sample ensemble member to file to check
  printf ("Writing out a sample member\n");
  sprintf (outfile, "Member12_%c", ::param);
  writefield ( ::outpath,
               outfile,
               &(ens[12]) );


  // Compute the mean and standard deviations for the large ensemble
  printf ("Computing the mean and standard deviation for the large ensemble\n");
  mean_stddev_calc ( ens,
                     ::nens,
                     &mean_full,
                     &std_full,
                     ::choosemembers93 );

  // Output full ensemble mean
  printf ("Outputting mean\n");
  sprintf (outfile, "FullMean_%c", ::param);
  writefield ( ::outpath,
               outfile,
               &mean_full );

  // Output full ensemble standard deviation
  printf ("Outputting standard deviation for full ensemble\n");
  sprintf (outfile, "FullStddev_%c", ::param);
  writefield ( ::outpath,
               outfile,
               &std_full );


  // Compute the mean and standard deviations for the reduced ensemble
  printf ("Computing the mean and standard deviation for the reduced ensemble\n");
  mean_stddev_calc ( ens,
                     ::nens,
                     &mean_red,
                     &std_red,
                     choosemembers );


  // Output reduced ensemble mean
  printf ("Outputting mean for the reduced ensemble\n");
  sprintf (outfile, "ReducedMean_%c_%02u", ::param, ::nens_red);
  writefield ( ::outpath,
               outfile,
               &mean_red );

  // Output reduced ensemble standard deviation
  printf ("Outputting standard deviation for the reduced ensemble\n");
  sprintf (outfile, "ReducedStd_%c_%02u", ::param, ::nens_red);
  writefield ( ::outpath,
               outfile,
               &std_red );



  // Go round each position in the domain
  for (x0=0; x0<(::Nx); x0++)
  { printf ("Long %u of %u\n", x0+1, ::Nx);
    for (y0=0; y0<(::Ny); y0++)

    { //printf ("  Lat %u of %u\n", y0+1, ::Ny);
      if (x0>(::radius) && x0<(::Nx-::radius) && y0>(::radius) && y0<(::Ny-::radius))

      { // Compute the correlation at this position
        x_low = x0 - ::radius;
        y_low = y0 - ::radius;
        x_hi  = x0 + ::radius;
        y_hi  = y0 + ::radius;

        //printf ("Computing the correlation field\n");
        init_field ( &corr );

        corr_calc ( ens,
                    ::nens,
                    &std_full,
                    &mean_red,
                    &corr,
                    x0,
                    y0,
                    x_low,
                    y_low,
                    x_hi,
                    y_hi,
                    choosemembers );

        // Output correlation field
/*
        printf ("Outputting correlation field\n");
        strcpy (outfile, "ExampleCorrelation");
        writesubfield ( ::outpath,
                        outfile,
                        &corr );
*/

        // Unravel this correlation field
        //printf ("Unravelling the correlation field\n");
        c = 0;
        for (x=0; x<(2*::radius+1); x++)
        { for (y=0; y<(2*::radius+1); y++)
          { corr1[c] = corr.value[x][y];
            c++;
          }
        }

        // Perform a fit to this correlation field

        // Starting values for the parameters
        apriori[0] = ens[0].xvals[x0];   // x0  (centre of ellipse)
        apriori[1] = ens[0].yvals[y0];   // y0  (centre of ellipse)

        //printf ("Performing least squares\n");
        LeastSq ( p,
                  apriori,
                  delta,
                  p_stddev,
                  y_stddev,
                  corr1,
                  corr.xvals,
                  corr.yvals,
                  ::maxiters,
                  critdiff );

        // Store the largest lengthscale (elements 2 and 3 are the lengthscales)
        if (p[2] > p[3])
        { lengths.value[x0][y0] = p[2];
        }
        else
        { lengths.value[x0][y0] = p[3];
        }
      }
      else
      { lengths.value[x0][y0] = 0.0;
      }

    }
  }

  // Output the lengtscales field
  printf ("Outputting the lengthscale file\n");
  sprintf (outfile, "Lengths_%c_%02d", ::param, nens_red);
  for (x0=0; x0<(::Nx); x0++)
  { lengths.xvals[x0] = ens[0].xvals[x0];
  }
  for (y0=0; y0<(::Ny); y0++)
  { lengths.yvals[y0] = ens[0].yvals[y0];
  }
  writefield ( ::outpath,
               outfile,
               &lengths );


  // Tidy up
  delete[] ens;

}



/* -------------------------------------------------------------------------------
   Extract key program arguments (model type and run type)
   ------------------------------------------------------------------------------- */
void ReadEns (
       char   path[],                  // Pathname
       char   param,                   // t or q
       int    t,                       // time index
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
  start[2] = 0;     count[2] = ::Ny   ;   // y
  start[1] = lev;   count[1] = 1;         // z
  start[0] = t;     count[0] = 1;         // t

  for (ensloop=0; ensloop<(::nens); ensloop++)
  { // Construct the filename
    sprintf (filename, "%s/qwq107.oper%02d.pp1.full.nc", path, ensloop);
    printf ("Reading file %s\n", filename);
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

// -------------------------------------------------------------------------------

void writefield ( char              path[],
                  char              file[],
                  struct field_type *field )
// Write field to a single file
{ //Declare local variables
  char         filename[::lenfilen];
  int          x, y;
  // netCDF-related variables
  size_t       Nx1, Ny1;
  int          ierr, ierr1, ierr2, ncid;
  int          dimidx, dimidy;
  int          varidx, varidy, varid;
  int          list1[1], list2[2];
  size_t       start[2], count[2];
  float        lineofdata[::Nx];


  sprintf (filename, "%s/%s.nc", path, file);
  printf ("Output filename : %s\n", filename);
  // Copy over dimension lengths to correct type
  Nx1 = ::Nx;
  Ny1 = ::Ny;

  // Create a new netCDF file
  ierr = nc_create (filename,
                    NC_CLOBBER,
                    &ncid);
  if (ierr != 0)
  { printf ("Error opening output file\n%s\nExiting\n", nc_strerror(ierr));
    exit(0);
  }

  // Define the dimensions
  ierr1 = nc_def_dim (ncid,
                      "x_distance",
                      Nx1,
                      &dimidx);
  ierr2 = nc_def_dim (ncid,
                      "y_distance",
                      Ny1,
                      &dimidy);
  if (ierr1 + ierr2 != 0)
  { printf ("%s\n%s\n", nc_strerror(ierr1), nc_strerror(ierr2));
    exit(0);
  }

  // Define the variables
  list1[0] = dimidx;
  ierr = nc_def_var (ncid,
                     "x_distance",
                     NC_FLOAT,
                     1,
                     list1,
                     &varidx);
  list1[0] = dimidy;
  ierr = nc_def_var (ncid,
                     "y_distance",
                     NC_FLOAT,
                     1,
                     list1,
                     &varidy);

  // x, y, z, pair directions correspond to 3, 2, 1, 0 respectively.
  list2[0] = dimidy;
  list2[1] = dimidx;

  ierr1 = nc_def_var (ncid,
                      "field",
                      NC_FLOAT,
                      2,
                      list2,
                      &varid);

  if (ierr1 != 0)
  { printf ("%s\n", nc_strerror(ierr1));
    exit(0);
  }

  // End define mode
  ierr = nc_enddef (ncid);
  if (ierr != 0)
  { printf ("Error changing mode\n%s\nExiting\n", nc_strerror(ierr));
    exit(0);
  }

  // Output data
  start[0] = 0;
  count[0] = Nx1;
  for (x=0; x<(::Nx); x++)
  { lineofdata[x] = float((*field).xvals[x]);
  }
  ierr1 = nc_put_var_float (ncid,
                            varidx,
                            lineofdata);
  count[0] = Ny1;
  for (y=0; y<(::Ny); y++)
  { lineofdata[y] = float((*field).yvals[y]);
  }
  ierr2 = nc_put_var_float (ncid,
                            varidy,
                            lineofdata);
  if (ierr1 + ierr2 != 0)
  { printf ("%s\n%s\n", nc_strerror(ierr1), nc_strerror(ierr2));
    exit(0);
  }

                    count[1] = 1;      // x
  start[0] = 0;     count[0] = Ny1;    // y

  for (x=0; x<(::Nx); x++)
  { start[1] = x;
    for (y=0; y<(::Ny); y++)
    { lineofdata[y] = float((*field).value[x][y]);
    }
    ierr1 = nc_put_vara_float (ncid,
                               varid,
                               start,
                               count,
                               lineofdata);
  }

  if (ierr1 != 0)
  { printf ("%s\n", nc_strerror(ierr1));
    exit(0);
  }

  //Close the netCDF file
  ierr = nc_close (ncid);
  if (ierr != 0)
  { printf ("Error closing write file\n%s\nExiting\n", nc_strerror(ierr));
    exit(0);
  }
}



// -------------------------------------------------------------------------------

void writesubfield ( char                  path[],
                     char                  file[],
                     struct sub_field_type *field )
// Write field to a single file
{ //Declare local variables
  char         filename[::lenfilen];
  int          x, y;
  // netCDF-related variables
  size_t       Nx1, Ny1;
  int          ierr, ierr1, ierr2, ncid;
  int          dimidx, dimidy;
  int          varidx, varidy, varid;
  int          list1[1], list2[2];
  size_t       start[2], count[2];
  float        lineofdata[2*::radius+1];


  sprintf (filename, "%s/%s.nc", path, file);
  printf ("Output filename : %s\n", filename);
  // Copy over dimension lengths to correct type
  Nx1 = 2*::radius+1;
  Ny1 = 2*::radius+1;

  // Create a new netCDF file
  ierr = nc_create (filename,
                    NC_CLOBBER,
                    &ncid);
  if (ierr != 0)
  { printf ("Error opening output file\n%s\nExiting\n", nc_strerror(ierr));
    exit(0);
  }

  // Define the dimensions
  ierr1 = nc_def_dim (ncid,
                      "x_distance",
                      Nx1,
                      &dimidx);
  ierr2 = nc_def_dim (ncid,
                      "y_distance",
                      Ny1,
                      &dimidy);
  if (ierr1 + ierr2 != 0)
  { printf ("%s\n%s\n", nc_strerror(ierr1), nc_strerror(ierr2));
    exit(0);
  }

  // Define the variables
  list1[0] = dimidx;
  ierr = nc_def_var (ncid,
                     "x_distance",
                     NC_FLOAT,
                     1,
                     list1,
                     &varidx);
  list1[0] = dimidy;
  ierr = nc_def_var (ncid,
                     "y_distance",
                     NC_FLOAT,
                     1,
                     list1,
                     &varidy);

  // x, y, z, pair directions correspond to 3, 2, 1, 0 respectively.
  list2[0] = dimidy;
  list2[1] = dimidx;

  ierr1 = nc_def_var (ncid,
                      "field",
                      NC_FLOAT,
                      2,
                      list2,
                      &varid);

  if (ierr1 != 0)
  { printf ("%s\n", nc_strerror(ierr1));
    exit(0);
  }

  // End define mode
  ierr = nc_enddef (ncid);
  if (ierr != 0)
  { printf ("Error changing mode\n%s\nExiting\n", nc_strerror(ierr));
    exit(0);
  }

  // Output data
  start[0] = 0;
  count[0] = Nx1;
  for (x=0; x<(2*::radius+1); x++)
  { lineofdata[x] = float((*field).xvals[x]);
  }
  ierr1 = nc_put_var_float (ncid,
                            varidx,
                            lineofdata);
  count[0] = Ny1;
  for (y=0; y<(2*::radius+1); y++)
  { lineofdata[y] = float((*field).yvals[y]);
  }
  ierr2 = nc_put_var_float (ncid,
                            varidy,
                            lineofdata);
  if (ierr1 + ierr2 != 0)
  { printf ("%s\n%s\n", nc_strerror(ierr1), nc_strerror(ierr2));
    exit(0);
  }

                    count[1] = 1;      // x
  start[0] = 0;     count[0] = Ny1;    // y

  for (x=0; x<(2*::radius+1); x++)
  { start[1] = x;
    for (y=0; y<(2*::radius+1); y++)
    { lineofdata[y] = float((*field).value[x][y]);
    }
    ierr1 = nc_put_vara_float (ncid,
                               varid,
                               start,
                               count,
                               lineofdata);
  }

  if (ierr1 != 0)
  { printf ("%s\n", nc_strerror(ierr1));
    exit(0);
  }

  //Close the netCDF file
  ierr = nc_close (ncid);
  if (ierr != 0)
  { printf ("Error closing write file\n%s\nExiting\n", nc_strerror(ierr));
    exit(0);
  }
}



// -------------------------------------------------------------------------------

void mean_stddev_calc ( struct field_type ens[],
                        int               nens,
                        struct field_type *mean,
                        struct field_type *std,
                        bool              choosemembers[] )
{ // Declare local variables
  int    x, y, mem;
  double nens_red, mn, var, dev;

  // Find out how many members are present in the reduced ensemble
  nens_red = 0.0;
  for (mem=0; mem<nens; mem++)
  { if (choosemembers[mem])
    { nens_red += 1.0;
    }
  }

  for (x=0; x<(::Nx); x++)
  { for (y=0; y<(::Ny); y++)
    { // Compute the mean
      mn = 0.0;
      for (mem=0; mem<nens; mem++)
      { if (choosemembers[mem])
        { mn += ens[mem].value[x][y];
        }
      }
      mn /= nens_red;
      (*mean).value[x][y] = mn;

      // Compute the variance
      var = 0.0;
      for (mem=0; mem<nens; mem++)
      { if (choosemembers[mem])
        { dev  = ens[mem].value[x][y] - mn;
          var += dev * dev;
        }
      }
      var /= (nens_red - 1.0);

      // Compute the standard deviation
      (*std).value[x][y] = sqrt(var);
    }
  }

  // Copy over the x and y values
  for (x=0; x<(::Nx); x++)
  { (*mean).xvals[x] = ens[0].xvals[x];
    (*std).xvals[x]  = ens[0].xvals[x];
  }
  for (y=0; y<(::Ny); y++)
  { (*mean).yvals[y] = ens[0].yvals[y];
    (*std).yvals[y]  = ens[0].yvals[y];
  }
}


// -------------------------------------------------------------------------------

void corr_calc ( struct field_type     ens[],
                 int                   nens,
                 struct field_type     *std,
                 struct field_type     *mean,
                 struct sub_field_type *corr,
                 int                   x0,     // The source point
                 int                   y0,
                 int                   x1,     // Lower point
                 int                   y1,
                 int                   x2,     // Upper point
                 int                   y2,
                 bool                  choosemembers[] )
{ // Declare local variables
  int    x, y, mem;
  double nens_red, cov, source_std;

  // Find out how many members are present in the reduced ensemble
  nens_red = 0.0;
  for (mem=0; mem<nens; mem++)
  { if (choosemembers[mem])
    { nens_red += 1.0;
    }
  }

  source_std = (*std).value[x0][y0];

  // Loop around the source point between the lower and upper points
  for (x=x1; x<=x2; x++)
  { for (y=y1; y<=y2; y++)
    { // Compute the covariance between (x0,y0) and (x,y)
      cov = 0.0;
      for (mem=0; mem<nens; mem++)
      { if (choosemembers[mem])
        { cov += (ens[mem].value[x0][y0] - (*mean).value[x0][y0]) *
                 (ens[mem].value[x][y] - (*mean).value[x][y]);
        }
      }
      cov /= nens_red;
      (*corr).value[x-x1][y-y1] = cov / (source_std * (*std).value[x][y]);
    }
  }

  // Copy over the x and y values
  for (x=x1; x<=x2; x++)
  { (*corr).xvals[x-x1] = (*mean).xvals[x];
  }
  for (y=y1; y<=y2; y++)
  { (*corr).yvals[y-y1] = (*mean).yvals[y];
  }

}

// -------------------------------------------------------------------------------

void init_field ( struct sub_field_type *field )

{ // Declare local variables
  int    x, y;
  for (x=0; x<(2*::radius+1); x++)
  { for (y=0; y<(2*::radius+1); y++)
    { (*field).value[x][y] = 0.0;
    }
  }
}



// -------------------------------------------------------------------------------

double exponential( double x,
                    double y,
                    double p[7] )  // x0, y0, Lxp, Lyp, C0, D0, rot

{ // Declare local variables
  double s, c;
  // Returns a gaussian function with the given parameters
  s = sin(p[6]);
  c = cos(p[6]);
  return p[4] * exp(-1.0 * fabs((x-p[0])*c + (y-p[1])*s) / p[2]) *
                exp(-1.0 * fabs((y-p[1])*c - (x-p[0])*s) / p[3]) + p[5];
}



// -------------------------------------------------------------------------------

void fn ( struct sub_field_type *field,
          double                p[7] )

{ // Declare local variables
  int    x, y;
  for (x=0; x<(2*::radius+1); x++)
  { for (y=0; y<(2*::radius+1); y++)
    { (*field).value[x][y] = exponential ( (*field).xvals[x],
                                           (*field).yvals[y],
                                           p );
    }
  }
}



// -------------------------------------------------------------------------------

void Fmatrix ( double xvals[2*::radius+1],
               double yvals[2*::radius+1],
               double p0[7],
               double delta[7],
               double F[(2*::radius+1)*(2*::radius+1)][7] )
// Compute the Jacobian
{ // Declare local variables
  int                   N, x, y, p1, p2, c;
  struct sub_field_type plus, minus;
  double                params[7];

  N = 7;

  // Copy over x and y values
  for (x=0; x<(2*::radius+1); x++)
  { plus.xvals[x]  = xvals[x];
    minus.xvals[x] = xvals[x];
  }
  for (y=0; y<(2*::radius+1); y++)
  { plus.yvals[y]  = yvals[y];
    minus.yvals[y] = yvals[y];
  }

  // Loop round each parameter to perturb
  for (p1=0; p1<N; p1++)
  { // Reset the parameters to their central values
    for (p2=0; p2<N; p2++)
    { params[p2] = p0[p2];
    }

    // Perturb this parameter negative
    params[p1] = p0[p1] - delta[p1];

    // Compute the function
    fn ( &minus,
         params );

    // Perturb this parameter positive
    params[p1] = p0[p1] + delta[p1];

    // Compute the function
    fn ( &plus,
         params );

    // Fill out the Jacobian
    c = 0;
    for (x=0; x<(2*::radius+1); x++)
    { for (y=0; y<(2*::radius+1); y++)
      { F[c][p1] = (plus.value[x][y] - minus.value[x][y]) / (2.0 * delta[p1]);
        c++;
      }
    }

  }
}



// -------------------------------------------------------------------------------

void FTF ( double F[2*::radius+1][7],
           double recip_y_stddev2,
           double product[7][7] )
// Compute the product F^T F
{ // Declare local variables
  int    c, p1, p2, len, N;
  double total;
  len = (2*::radius+1) * (2*::radius+1);
  N   = 7;

  // Deal with the upper-triangular part
  for (p1=0; p1<N; p1++)
  { for (p2=p1; p2<N; p2++)
    { total = 0.0;
      for (c=0; c<len; c++)
      { total += F[c][p1] * F[c][p2];
      }
      product[p1][p2] = recip_y_stddev2 * total;
    }
  }

  // Use symmetry to deal with the remaining part
  for (p1=1; p1<N; p1++)
  { for (p2=0; p2<p1; p2++)
    { product[p1][p2] = product[p2][p1];
    }
  }

}



// -------------------------------------------------------------------------------

void GaussEl ( double A[7][7],       // Matrix (modified on output)
               double x[7],          // Solution (output)
               double y[7] )         // RHS (modified on output)
// Gaussian elimination for a full square matrix
{ // Declare local variables
  int    N=7;
  int    i, j, jp;
  double recip_Aii, Aii, Aji, partialsum;

  // First go through the matrix
  for (i=0; i<(N-1); i++)
  { recip_Aii = 1.0 / A[i][i];
    // Modify later rows accordingly
    for (j=i+1; j<N; j++)
    { Aji = A[j][i];
      // Modify RHS
      y[j] -= Aji * y[i] * recip_Aii;
      // Modify matrix
      for (jp=i+1; jp<N; jp++)
      { A[j][jp] -= Aji * A[i][jp] * recip_Aii;
      }
    }
  }

  // Now solve the matrix backwards
  for (i=N-1; i>=0; i--)
  { partialsum = 0.0;
    if (i<(N-1))
    { // Compute partial sum
      for (j=i+1; j<N; j++)
      { partialsum += A[i][j] * x[j];
      }
    }
    Aii       = A[i][i];
    recip_Aii = 1.0 / Aii;
    x[i]        = (y[i] - partialsum) * recip_Aii;
  }
}

// -------------------------------------------------------------------------------

void LeastSq ( double p[7],                // Parameters to optimise
               double apriori[7],          // a-priori parameter values
               double delta[7],            // Small perts to parameters (for finite differencing)
               double p_stddev[7],         // sqroot of diagonal of P
               double y_stddev,            // sqroot of diagonal of R
               double y[(2*::radius+1)*(2*::radius+1)],     // Data to fit to
               double xvals[2*::radius+1], // x-values (co-ordinates)
               double yvals[2*::radius+1], // y-values (co-ordinates)
               int    maxiters,            // Max number of iterations
               double critdiff )           // Stopping criteria

// Parameter fitting
{ // Declare local variables
  int                   iter, c, len, N, p1, p2, xx, yy;
  struct sub_field_type field;
  double                deltap_mag, total;
  double                p_hat[7], deltap[7], deltapb[7];
  double                Hessian[7][7];
  double                diff[(2*::radius+1)*(2*::radius+1)];
  double                F[(2*::radius+1)*(2*::radius+1)][7];
  double                recip_y_stddev2, recip_p_stddev2[7];
  // Temporary stuff (for testing)
  //char                  filename[::lenfilen];
  //FILE*                 op = NULL;

  //sprintf (filename, "%s/Fit.dat", ::outpath);
  //printf ("Output filename : %s\n", filename);
  //op = fopen (filename, "w");
  //fprintf (op, "#iter   x0        y0        Lxp        Lyp        C0        D0        rot\n");


  len  = (2*::radius+1) * (2*::radius+1);
  iter = 0;
  N    = 7;

  // Preprocess covariance stats for efficiency
  recip_y_stddev2 = 1.0 / (y_stddev * y_stddev);
  for (p1=0; p1<N; p1++)
  { recip_p_stddev2[p1] = 1.0 / (p_stddev[p1] * p_stddev[p1]);
  }

  // Copy over x and y values
  for (xx=0; xx<(2*::radius+1); xx++)
  { field.xvals[xx]  = xvals[xx];
  }
  for (yy=0; yy<(2*::radius+1); yy++)
  { field.yvals[yy]  = yvals[yy];
  }

  // Start at the a-priori value
  //printf ("      x0     y0     Lxp     Lyp     C0     D0     rot\n");
  //printf ("pb:  ");
  //fprintf (op, "%u  ", iter);
  for (p1=0; p1<N; p1++)
  { p[p1]       = apriori[p1];
    //printf ("%f  ", p[p1]);
    //fprintf (op, "%f  ", p[p1]);
  }
  //printf ("\n");
  //fprintf (op, "\n");


  do
  { //printf ("Starting iteration %u\n", iter);
    //printf ("  Compting the model data\n");

    // Evaluate deltapb
    for (p1=0; p1<N; p1++)
    { deltapb[p1] = apriori[p1] - p[p1];
    }

    // Compute the model's version of the data
    fn ( &field,
         p );

    // Output the field produced by this set of parameters
    //printf ("  Outputting the model data\n");
    //sprintf (filename, "Iteration%02d", iter);
    //writefield ( ::outpath,
    //             filename,
    //             &field );

    // Compute the scaled difference with the actual data
    //printf ("  Computing the scaled difference\n");
    c = 0;
    for (xx=0; xx<(2*::radius+1); xx++)
    { for (yy=0; yy<(2*::radius+1); yy++)
      { diff[c] = recip_y_stddev2 * (y[c] - field.value[xx][yy]);
        c++;
      }
    }

    // Compute F
    //printf ("  Computing the Jacobian, F\n");
    Fmatrix ( xvals,
              yvals,
              p,
              delta,
              F );

    // Act with F^T
    //printf ("  Acting with F^T\n");
    for (p1=0; p1<N; p1++)
    { total = 0.0;
      for (c=0; c<len; c++)
      { total += F[c][p1] * diff[c];
      }
      p_hat[p1] = total;
    }

    // Add on term associated with prior
    for (p1=0; p1<N; p1++)
    { p_hat[p1] += recip_p_stddev2[p1] * deltapb[p1];
    }

    // Construct F^T R^-1 F
    //printf ("  Constructing F^T R-1 F\n");
    FTF ( F,
          recip_y_stddev2,
          Hessian );


    // Add on a-priori bit to the Hessian
    for (p1=0; p1<N; p1++)
    { Hessian[p1][p1] += recip_p_stddev2[p1];
    }

    // Print out the Hessian
    //printf ("  ---- Hessian ----------------\n");
    //for (p1=0; p1<N; p1++)
    //{ printf ("  %u:  ", p1);
    //  for (p2=0; p2<N; p2++)
    //  { printf ("%f ", Hessian[p1][p2]);
    //  }
    //  printf ("\n");
    //}
    //printf ("  -----------------------------\n");

    // Solve the linear system for the deltap vector
    //printf ("  Solving the linear system\n");
    GaussEl ( Hessian,
              deltap,
              p_hat );

    // Modify p
    //printf ("  Modifying p\n");
    //printf ("      x0     y0     Lxp     Lyp     C0     D0     rot\n");
    //printf ("  p:  ");
    //fprintf (op, "%u  ", iter);
    for (p1=0; p1<N; p1++)
    { p[p1] += deltap[p1];
      //printf ("%f  ", p[p1]);
      //fprintf (op, "%f  ", p[p1]);
    }
    //printf ("\n");
    //fprintf (op, "\n");

    // Compute the size of deltap
    deltap_mag = 0.0;
    for (p1=0; p1<N; p1++)
    { deltap_mag += recip_p_stddev2[p1] * deltap[p1] * deltap[p1];
    }
    deltap_mag = sqrt(deltap_mag);
    //printf ("  Size of deltap at the end of iteration %u: %f\n", iter, deltap_mag);

    iter++;
  }
  while ((iter < maxiters) && (deltap_mag > critdiff));

  //fclose (op);

}
