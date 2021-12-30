#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"

using namespace CLHEP;

static const unsigned flatentries = 2;
static const double minenergy = 1.0 * eV;
static const double maxenergy = 8.0 * eV;

G4OpticalSurface *
MakeS_TyvekCrystal()
{
  const unsigned num = flatentries;
  double Ephoton[num] = {minenergy, maxenergy};
  double Reflectivity[num] = {0.979, 0.979};

  //////////////////////////////////
  // Realistic Crystal-Tyvek surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("TyvekOpSurface");
  surface->SetType(dielectric_LUT);
  surface->SetModel(LUT);
  surface->SetFinish(polishedtyvekair);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();

  table->AddProperty("REFLECTIVITY", Ephoton, Reflectivity, num);

  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_ESR()
{
  const unsigned num = flatentries;
  double Ephoton[num] = {minenergy, maxenergy};
  double Reflectivity[num] = {0.985, 0.985};
  // source: https://www.osti.gov/servlets/purl/1184400
  //////////////////////////////////
  // ESR surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("ESROpSurface");
  surface->SetType(dielectric_LUT);
  surface->SetModel(LUT);
  surface->SetFinish(polishedvm2000air);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();

  table->AddProperty("REFLECTIVITY", Ephoton, Reflectivity, num);

  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_IdealTyvekCrystal()
{
  //////////////////////////////////
  // Ideal Crystal-Tyvek surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("IdealTyvekOpSurface");
  surface->SetType(dielectric_LUT);
  surface->SetModel(LUT);
  surface->SetFinish(polishedtyvekair);

  return surface;
}

G4OpticalSurface *
MakeS_Polished()
{

  static const unsigned nentries = flatentries;
  static double phoE[nentries] = {minenergy, maxenergy};
  static double specularlobe[nentries] = {1.0, 1.0};

  //////////////////////////////////
  // Realistic polished surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("PolishedOpSurface");
  surface->SetType(dielectric_dielectric);
  surface->SetModel(unified);
  surface->SetFinish(ground);
  // necessary even for polished surfaces to enable UNIFIED code
  surface->SetSigmaAlpha(1.3 * degree); // Janecek2010 (1.3 * degree)

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();

  table->AddProperty(
      "SPECULARLOBECONSTANT", phoE, specularlobe, nentries);
  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_IdealPolished()
{
  //////////////////////////////////
  // Ideal polished surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("IdealOpSurface");
  surface->SetType(dielectric_dielectric);
  surface->SetModel(glisur);
  surface->SetFinish(polished);

  return surface;
}

G4OpticalSurface *
MakeS_Mirror()
{

  const unsigned nentries = flatentries;
  double phoE[nentries] = {minenergy, maxenergy};
  double reflectivity[nentries] = {0.9, 0.9};
  //////////////////////////////////
  // Mirror surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("MirrorOpSurface");
  surface->SetType(dielectric_metal);
  surface->SetFinish(polished);
  surface->SetModel(unified);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("REFLECTIVITY", phoE, reflectivity, nentries);
  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_IdealMirror()
{
  const unsigned nentries = flatentries;
  double phoE[nentries] = {minenergy, maxenergy};
  double reflectivity[nentries] = {1.0, 1.0};
  //////////////////////////////////
  // Ideal mirror surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("MirrorOpSurface");
  surface->SetType(dielectric_metal);
  surface->SetFinish(polished);
  surface->SetModel(unified);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("REFLECTIVITY", phoE, reflectivity, nentries);
  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_IdealWhiteSurface()
{
  //////////////////////////////////
  // Ideal mirror surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("WhiteOpSurface");
  surface->SetType(dielectric_metal);
  surface->SetFinish(ground);
  surface->SetModel(unified);

  double phoE[flatentries] = {minenergy, maxenergy};
  double reflectivity[flatentries] = {0.5, 0.5};

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("REFLECTIVITY", phoE, reflectivity, flatentries);
  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_Absorbing()
{
  const unsigned nentries = flatentries;
  double phoE[nentries] = {minenergy, maxenergy};
  double reflectivity[nentries] = {0.0, 0.0};
  //////////////////////////////////
  // Absorbing surface
  //////////////////////////////////
  G4OpticalSurface *surface = new G4OpticalSurface("AbsorbingOpSurface");
  surface->SetType(dielectric_dielectric);
  surface->SetFinish(groundfrontpainted);
  surface->SetModel(unified);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("REFLECTIVITY", phoE, reflectivity, nentries);
  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_wrap()
{
  //////////////////////////////////
  // wrap surface
  //////////////////////////////////
  G4OpticalSurface *scintWrap = new G4OpticalSurface("ScintWrap");
  scintWrap->SetType(dielectric_metal);
  scintWrap->SetFinish(ground);
  scintWrap->SetModel(unified);
  double pp[2] = {0.1 * eV, 6 * eV};
  double reflectivity[2] = {1, 1};
  G4MaterialPropertiesTable *scintWrapProperty = new G4MaterialPropertiesTable();
  scintWrapProperty->AddProperty("REFLECTIVITY", pp, reflectivity, 2);
  scintWrap->SetMaterialPropertiesTable(scintWrapProperty);
  return scintWrap;
}

G4OpticalSurface *MakeS_PMT()
{

  const unsigned nentries = 21;
  double phoE[nentries] = {4.591666802 * eV, 4.249963562 * eV, 3.970289068 * eV, 3.545647077 * eV,
                           3.157143617 * eV, 2.8391158 * eV, 2.630406004 * eV, 2.436700346 * eV, 2.337427879 * eV, 2.246008967 * eV,
                           2.193827957 * eV, 2.102995234 * eV, 2.041160762 * eV, 1.99177784 * eV, 1.910653049 * eV, 1.872420367 * eV,
                           1.830538262 * eV, 1.812740321 * eV, 1.78790683 * eV, 1.768575074 * eV, 1.761389217 * eV};
  double efficiency[nentries] = {13.033, 20.201, 23.97, 25.352, 25.317, 22.981, 18.968, 15.068,
                                 9.168, 6.3757, 5.1664, 3.3926, 2.2282, 1.3305, 0.40721, 0.18973, 
                                 0.073042, 0.043624, 0.021119, 0.012856, 0.010031};
  for (int dd = 0; dd < 21; dd++)
    efficiency[dd] /= 100.;

  const unsigned ref_ent = 2;
  double phoE2[ref_ent] = {1 * eV, 6 * eV};
  double refraction[ref_ent] = {1.473, 1.473};
  double QE[ref_ent] = {1.0, 1.0};

  G4OpticalSurface *surface = new G4OpticalSurface("PMT_Surface");
  surface->SetType(dielectric_metal);
  surface->SetModel(unified);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();
  table->AddProperty("EFFICIENCY", phoE2, QE, ref_ent);

  double reflectivity[ref_ent] = {0.0, 0.0};

  table->AddProperty("REFLECTIVITY", phoE2, reflectivity, ref_ent);
  surface->SetMaterialPropertiesTable(table);

  return surface;
}

G4OpticalSurface *
MakeS_PCBSurface()
{
  // PCB is a flat gray surface for now
  double phoE2[flatentries] = {minenergy, maxenergy};
  double reflectivity[flatentries] = {0.5, 0.5};

  G4OpticalSurface *surface = new G4OpticalSurface("PCB_Surface");
  surface->SetType(dielectric_metal);
  surface->SetFinish(ground);
  surface->SetModel(unified);

  G4MaterialPropertiesTable *table = new G4MaterialPropertiesTable();

  table->AddProperty("REFLECTIVITY", phoE2, reflectivity, flatentries);
  surface->SetMaterialPropertiesTable(table);
  return surface;
}
