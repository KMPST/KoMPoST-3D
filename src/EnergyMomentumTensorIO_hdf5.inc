#include <cstdlib>
#include <cstring>
#include <gsl/gsl_math.h>
#include <hdf5.h>
#include <iostream>
#include <string>
#include <vector>

void write_hdf5_slice(hid_t dataset, float *buffer, int rank,
                      hsize_t size[], hsize_t offset[]) {
  hid_t file_space=H5Dget_space(dataset);
  hid_t memory_space=H5Screate_simple(rank,size,NULL);
  H5Sselect_hyperslab(file_space,H5S_SELECT_SET,offset,NULL,size,NULL);
  if (H5Dwrite(dataset,H5T_NATIVE_FLOAT,memory_space,file_space,
               H5P_DEFAULT,buffer)<0) exit(EXIT_FAILURE);
  H5Sclose(memory_space);
  H5Sclose(file_space);
}

void read_hdf5_slice(hid_t dataset, float *buffer, int rank,
                     hsize_t size[], hsize_t offset[]) {
  hid_t file_space=H5Dget_space(dataset);
  hid_t memory_space=H5Screate_simple(rank,size,NULL);
  H5Sselect_hyperslab(file_space,H5S_SELECT_SET,offset,NULL,size,NULL);
  if (H5Dread(dataset,H5T_NATIVE_FLOAT,memory_space,file_space,
              H5P_DEFAULT,buffer)<0) exit(EXIT_FAILURE);
  H5Sclose(memory_space);
  H5Sclose(file_space);
}

void EnergyMomentumTensorMapLoadMcDipperHDF5(
    EnergyMomentumTensorMap *TIn, ChargeCurrentMap *JIn,
    std::string filename) {
  std::cerr << "#GETTING McDIPPER ENERGY AND CHARGES FROM HDF5 FILE "
            << filename << std::endl;
  hid_t file=H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  if (file<0) exit(EXIT_FAILURE);
  hid_t dataset=H5Dopen2(file,"charges",H5P_DEFAULT);
  if (dataset<0) exit(EXIT_FAILURE);
  hid_t data_space=H5Dget_space(dataset);
  hsize_t dimensions[2];
  const int rank=H5Sget_simple_extent_ndims(data_space);
  H5Sget_simple_extent_dims(data_space,dimensions,NULL);
  const hsize_t number_of_cells=(hsize_t)EventInput::Neta
      *EventInput::Ns*EventInput::Ns;
  if (rank!=2 || dimensions[0]!=number_of_cells || dimensions[1]<5) {
    std::cerr << "Error: incompatible McDIPPER /charges shape" << std::endl;
    exit(EXIT_FAILURE);
  }
  H5Sclose(data_space);
  const int cells_per_slice=EventInput::Ns*EventInput::Ns;
  std::vector<float> density_slice(cells_per_slice*5);
   const double unit_scale=M_HBARC*M_HBARC*M_HBARC;
  const double small_num=1000*GSL_DBL_MIN;
  double charge_scale=0.0;
  if (JIn!=NULL) charge_scale=unit_scale/JIn->tau;
  for (int ieta=0;ieta<EventInput::Neta;ieta++) {
    hsize_t size[2]={(hsize_t)cells_per_slice,5};
    hsize_t offset[2]={(hsize_t)ieta*cells_per_slice,
                       dimensions[1]-5};
    read_hdf5_slice(dataset,density_slice.data(),2,size,offset);
    for (int ix=0;ix<EventInput::Ns;ix++)
      for (int iy=0;iy<EventInput::Ns;iy++) {
        const int index=(ix*EventInput::Ns+iy)*5;
        const double tau_ed=KoMPoSTInputs::McDipperKFactor
                            *density_slice[index]+density_slice[index+1];
        const double energy_density=tau_ed/TIn->tau;
        const double value=unit_scale*energy_density;
        if (energy_density<small_num)
          TIn->Set(ix,iy,ieta,10*small_num,5*small_num,5*small_num,
                   0.0,0.0,0.0,0.0,0.0,0.0,0.0);
        else
          TIn->Set(ix,iy,ieta,value,value/2.0,value/2.0,0.0,
                   0.0,0.0,0.0,0.0,0.0,0.0);
        if (JIn!=NULL) {
          JIn->Set(ChargeIndex::Up,ix,iy,ieta,
                   charge_scale*density_slice[index+2],0.0,0.0,0.0);
          JIn->Set(ChargeIndex::Down,ix,iy,ieta,
                   charge_scale*density_slice[index+3],0.0,0.0,0.0);
          JIn->Set(ChargeIndex::Strange,ix,iy,ieta,
                   charge_scale*density_slice[index+4],0.0,0.0,0.0);
        }
      }
  }
  H5Dclose(dataset);
  H5Fclose(file);
}

hid_t create_hdf5_dataset(hid_t file, const char *name, int rank,
                          hsize_t dimensions[]) {
  hid_t data_space=H5Screate_simple(rank,dimensions,NULL);
  hid_t dataset=H5Dcreate2(file,name,H5T_IEEE_F32LE,data_space,
                           H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  H5Sclose(data_space);
  if (dataset<0) exit(EXIT_FAILURE);
  return dataset;
}

void write_hdf5_double_attribute(hid_t file, const char *name,
                                 double value) {
  hid_t data_space=H5Screate(H5S_SCALAR);
  hid_t attribute=H5Acreate2(file,name,H5T_IEEE_F64LE,data_space,
                             H5P_DEFAULT,H5P_DEFAULT);
  if (attribute<0 ||
      H5Awrite(attribute,H5T_NATIVE_DOUBLE,&value)<0) exit(EXIT_FAILURE);
  H5Aclose(attribute);
  H5Sclose(data_space);
}

void write_hdf5_int_attribute(hid_t file, const char *name, int value) {
  hid_t data_space=H5Screate(H5S_SCALAR);
  hid_t attribute=H5Acreate2(file,name,H5T_STD_I32LE,data_space,
                             H5P_DEFAULT,H5P_DEFAULT);
  if (attribute<0 ||
      H5Awrite(attribute,H5T_NATIVE_INT,&value)<0) exit(EXIT_FAILURE);
  H5Aclose(attribute);
  H5Sclose(data_space);
}

void write_hdf5_string_attribute(hid_t object, const char *name,
                                 const char *value) {
  hid_t data_type=H5Tcopy(H5T_C_S1);
  H5Tset_size(data_type,strlen(value)+1);
  hid_t data_space=H5Screate(H5S_SCALAR);
  hid_t attribute=H5Acreate2(object,name,data_type,data_space,
                             H5P_DEFAULT,H5P_DEFAULT);
  if (attribute<0 || H5Awrite(attribute,data_type,value)<0)
    exit(EXIT_FAILURE);
  H5Aclose(attribute);
  H5Sclose(data_space);
  H5Tclose(data_type);
}

hid_t create_hdf5_output_file(std::string filename, double tau,
                              double dx, double dy, double deta) {
  hid_t file=H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,
                       H5P_DEFAULT,H5P_DEFAULT);
  if (file<0) exit(EXIT_FAILURE);
  write_hdf5_int_attribute(file,"format_version",1);
  write_hdf5_double_attribute(file,"tau_fm",tau);
  write_hdf5_double_attribute(file,"dx_fm",dx);
  write_hdf5_double_attribute(file,"dy_fm",dy);
  write_hdf5_double_attribute(file,"deta",deta);
  return file;
}

void create_hdf5_output_datasets(hid_t file, int Neta, int Nx, int Ny,
                                 int Ncharge, hid_t &hydro_dataset,
                                 hid_t &charge_dataset) {
  hsize_t hydro_size[4]={(hsize_t)Neta,(hsize_t)Nx,(hsize_t)Ny,16};
  hydro_dataset=create_hdf5_dataset(file,"hydro",4,hydro_size);
  write_hdf5_string_attribute(hydro_dataset,"axes","eta,x,y,field");
  write_hdf5_string_attribute(
      hydro_dataset,"fields",
      "energy_density,u_tau,u_x,u_y,u_eta,"
      "pi_tautau,pi_taux,pi_tauy,pi_taueta,pi_xx,"
      "pi_xy,pi_xeta,pi_yy,pi_yeta,pi_etaeta,P_iso");
  charge_dataset=-1;
  if (Ncharge>0) {
    hsize_t charge_size[5]={(hsize_t)Neta,(hsize_t)Nx,(hsize_t)Ny,
                            (hsize_t)Ncharge,5};
    charge_dataset=create_hdf5_dataset(file,"charge",5,charge_size);
    write_hdf5_string_attribute(
        charge_dataset,"axes","eta,x,y,charge,field");
    write_hdf5_string_attribute(charge_dataset,"charges","u,d,s");
    write_hdf5_string_attribute(
        charge_dataset,"fields","n,V_tau,V_x,V_y,V_eta");
  }
}

void close_hdf5_output(hid_t file, hid_t hydro_dataset,
                       hid_t charge_dataset) {
  H5Dclose(hydro_dataset);
  if (charge_dataset>=0) H5Dclose(charge_dataset);
  H5Fclose(file);
}
