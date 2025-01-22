// NOTE: The RNTuple classes are experimental at this point.
// Functionality, interface, and data format is still subject to changes.
// Do not use for real data!
 
#include <ROOT/RNTupleDS.hxx>
#include <ROOT/RNTupleImporter.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RPageStorageFile.hxx>
 
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
 
// Import classes from experimental namespace for the time being.
// follows https://root.cern/doc/master/ntpl008__import_8C.html
// tested in 6.35.01 
using RNTuple = ROOT::RNTuple;
using RNTupleImporter = ROOT::Experimental::RNTupleImporter;
using RNTupleReader = ROOT::Experimental::RNTupleReader;
 
constexpr char const *kTreeName = "Events";

void RNtuples(const char* kTreeFileName, const char* kNTupleFileName)
{

   // RNTupleImporter appends keys to the output file; make sure a second run of the tutorial does not fail
   // with `Key 'Events' already exists in file ntpl008_import.root` by removing the output file.
   gSystem->Unlink(kNTupleFileName);
 
   // Use multiple threads to compress RNTuple data.
   ROOT::EnableImplicitMT();
 
   // Create a new RNTupleImporter object.
   auto importer = RNTupleImporter::Create(kTreeFileName, kTreeName, kNTupleFileName);

   // By default, the RNTuple is compressed with zstd
   // to compress the imported RNTuple using lz4 (with compression level 4) instead 
   auto writeOptions = importer->GetWriteOptions(); 
   writeOptions.SetCompression(404);
   importer->SetWriteOptions(writeOptions);
   
   // Begin importing.
   importer->Import();
 
   // Inspect the schema of the written RNTuple.
   auto file = std::unique_ptr<TFile>(TFile::Open(kNTupleFileName));
   if (!file || file->IsZombie()) {
      std::cerr << "cannot open " << kNTupleFileName << std::endl;
      return;
   }
   auto ntpl = std::unique_ptr<RNTuple>(file->Get<RNTuple>("Events"));
   auto reader = RNTupleReader::Open(*ntpl);
   reader->PrintInfo();

}
