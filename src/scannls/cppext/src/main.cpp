#include "rescuer.h"
// solution for segfault in future
// https://github.com/corwinjoy/Complete-Striped-Smith-Waterman-Library/pull/1/commits/10de7dbb01d04f9859571d79f3ecc633ebeff99e

int main(int argc, char* argv[]) {
  std::cout << argc << argv[0] << std::endl;

  std::string file{"/panfs/home/yang4414/li002252/project/scan_data/ONT_PC3_chr1_10.bam"};

  int mapq_threshold = 15;
  int min_softclip_length = 5;
  int min_mismatch_count = 3;
  double min_align_ratio = 0.8;

  std::vector<std::string> current_names{
      "b34721e7-34f1-451d-b978-9c0c092092e3", "44f31eb5-e224-4e7d-85d9-f59d16915df0",
      "6f25210a-04e0-4969-9402-13b088db558f", "2372e92c-6ac9-4156-884d-9e1c767575e2",
      "80acdb57-c0fa-4fff-8cac-4395c89110cf", "87d34704-fe39-492f-be2e-e24963d87c77",
      "c67f7ed5-1585-40de-bb93-da07f5b7f6d1", "8bc5ed50-339d-4d56-9327-c475b71314dd",
      "ccc6ec1f-46fa-412a-9943-160243728417", "9578af32-d82d-45dc-92c1-3ec0aa265ad1"};

  std::vector<std::string> names_list{
      "62b40e84-93ea-4b98-92af-a48d962a1366", "16004482-b8a3-4c33-84b9-5b17b4bb6c2f",
      "0538cbf6-01e4-4017-afd0-29b5731f1fd0", "24d89428-1437-4cc0-98d0-049d5982a3ab",
      "59505063-87ce-445f-b6d9-0c612502db2e", "1b6a6723-0809-4d40-b5d3-37704f5923df",
      "73686bbc-b255-4978-acd3-d7809c6ec252", "d31c3d90-1b00-4297-ac8b-284f64511a41",
      "202bcbea-f41a-4bbd-890a-ba9af8be20f1", "6da140db-1c0a-4dd9-a4d6-05b7b01460f9",
      "ccc6ec1f-46fa-412a-9943-160243728417", "cf26b46a-4c9e-46b1-8637-5eed458f566b",
      "78c22efd-48b3-4a55-b740-fc48f2f678c5", "bf9c4f15-1a91-4c38-91cd-298ed638b328",
      "c67f7ed5-1585-40de-bb93-da07f5b7f6d1", "86eaaf18-9493-48e6-b913-94361798e89f",
      "5e797ef4-d867-4d16-8707-3669ab06f95f", "11ddcdc6-0783-4a90-a5c3-73c82bd83740",
      "12ee9839-955f-46e5-99e5-07ba90aba8b4", "8999614e-885f-4c05-8e5f-ace1a86bd11e",
      "b85026de-416e-4667-854a-bfacd55c1956", "879ff446-4f88-4239-a547-e9b2e5ddb8c8",
      "3f76ad5b-8828-4dcf-b778-7b990f101690", "d2159eb7-823a-4903-90ed-10a72e6cd482",
      "cfc657f9-d8aa-47ae-9a7e-d2462073ea06", "29cf0325-771d-41b2-96e3-f02923667cc1",
      "3d4392c4-e812-45ad-8256-8808b2f2ddbe", "15c10cc4-038b-4e6b-b920-0a49d9f67251",
      "87d34704-fe39-492f-be2e-e24963d87c77", "8b45472d-e3bc-4a53-a41a-6e3a58e7c5e5",
      "d41e7e4f-1405-485e-8084-33ecbd4b833f", "b34721e7-34f1-451d-b978-9c0c092092e3",
      "bdf1c931-bbc1-461c-bfed-f11070fd506c", "cd1114c8-856b-46a7-9825-d6436cf878cb",
      "ce4a208f-5d98-4858-9874-dbd86667ba58", "5b7c228a-ec99-420d-b1c4-37db9b9ff850",
      "d2b2ae45-b9f2-4d50-8b41-b57f46d3ba93", "37512109-6383-44f1-9a8d-b23ec8f25e5b",
      "8d511a2f-e112-46f0-982c-797b4dbf67cc", "03e0716a-851c-48d2-b5c3-5dabe03414b4",
      "f28c5e7e-e12e-4f27-bab0-5cae2ba92d66", "9fcdef6a-6f41-4256-9ffd-fc28f95f00a7",
      "0d3f349b-f3b3-4cbf-ada2-d156109432ca", "c7d799bf-89ef-4aec-be18-f8640d0b1538",
      "0ddbf001-2f1d-4b9c-87cb-380fe7df73ca", "83a70225-d1dd-4c8d-93ba-477edc376dbc",
      "b7734ef0-75a5-45b4-b5b3-b2e7b2e4a3f6", "b22f29d2-86ee-4f60-8aac-fb3ab72dd13b",
      "932faec1-b378-435c-9ff3-12598b8e1e8e", "df4d35df-5a6b-4eb8-b695-eeeea264fc75",
      "b5278934-e936-4642-963b-0a461c8696d1", "08a3cf2e-e8a3-4dec-bd6b-4ae11ca90bc6",
      "11258c25-3dce-4d6f-9ebb-6050b95e5914", "0a41b2a6-8944-4568-9dbe-1161c2e3d8d2",
      "08e2a66e-975e-47c3-ac22-9954f5604021", "8fcfe228-517c-475b-a9a2-b29b01ae3832",
      "3c2524fa-98ef-476c-8182-6a5e733a4a30", "4385c532-13b1-4a5e-a1bd-24e447406fad",
      "296b12e5-777b-4ba1-9a62-336b27a79f37", "56dbbde5-a5d7-4f2b-bbaf-eeff2523dd2f",
      "112471a4-cbe8-4b07-835f-24a16dc13fe7", "c1cec480-d592-4dfd-8b1e-bbe0b1cb8072",
      "7a593f6c-3f8d-451d-bcfe-bd075ab031b0", "6f25210a-04e0-4969-9402-13b088db558f",
      "0ec91d98-603a-4084-87f9-5a720774f935", "9506cef1-d8e5-49cc-a0bd-44d6c2f1d07a",
      "9c2d608f-83c7-40a3-bbfa-7b58d0347c6d", "502fdd4a-f29c-43cc-8c7c-a51a5f4e0921",
      "e93c6403-a4c6-4992-9461-a8784f93c4a0", "ed4e407b-27b5-4ce0-8a2d-24adb3b776b8",
      "b8a005a0-44c2-4401-8961-e7387af63668", "80acdb57-c0fa-4fff-8cac-4395c89110cf",
      "8a705dea-a150-44fc-9d84-450a1354b947", "58e7f9f3-082b-44e4-8d85-87ba1c7158f0",
      "2f6c3cea-fbe8-43e1-8ec7-1fe1fee197e7", "a7e0719d-9f14-43c4-bc61-d37ed0991894",
      "3653e553-465e-4a21-b54b-b70dc655e4d6", "57a792e1-51a2-4109-81c0-656b2d0dea3b",
      "134f0c51-1bec-47d4-beff-ae480cf67067", "44f31eb5-e224-4e7d-85d9-f59d16915df0",
      "71ce86e2-583e-436d-a886-5911c2759972", "504ad005-0e64-4006-8a3c-441a458a0a47",
      "2efece49-2894-452b-be77-7d20466e279d", "f4196290-4068-4cbf-a89d-33c533a0d006",
      "53a19d3b-7a73-40ed-aadd-24f11b075d32", "2372e92c-6ac9-4156-884d-9e1c767575e2",
      "9578af32-d82d-45dc-92c1-3ec0aa265ad1", "e876055f-b702-4e53-b758-42a8724a7a2d",
      "624eda85-0ea7-4641-967d-cba874cb17b6", "779b5ba4-bd44-43a5-8699-b358d27b13f4",
      "4fdebab7-26c9-4693-8e6d-87e3c22c9322", "8bc5ed50-339d-4d56-9327-c475b71314dd",
      "5856831d-a0be-43a6-a54e-c2ae3cad6c52"};

  rescuer::Rescuer p_rescuer{file.c_str(),       mapq_threshold,  min_softclip_length,
                             min_mismatch_count, min_align_ratio, 10};

  int sr{-1};

  //    chrom = 'chr1' start = 199782025 mode2 = 2

  long start{199782025};
  std::string strand{"+"};
  sr = p_rescuer.calculate_sr("chr1", start, start, 2, strand, current_names,
                              names_list);  // chr17:7708250-7708250

  //  std::string query{
  //      "AGAAAGACTTTTCACAGAACAGACTATTACAGTATACTCGGGACACCATCTCTTACGATATCTTTAAAACCTAGGTGTCTGATTTCATGCTC"
  //      "TCCCTTAAAAAGTGTTCCTCTACCAACTATGAACAGGAATCATAGTCCTGTTAGAGAT"};
  //  std::string reference{
  //      "AGAAAGATTTCACAGAACAGACAATTACAGTACTCGGGGGGAACACCACCATCTTCTTACGATATTTAAAACCTAGGTGTCTTGATTTTCAT"
  //      "GCTACCTCCTTAAAAAGTGTCCTCTCTACCAACTATGAACAGAGTTCATAGTCCTGTT"};
  //
  //  std::cout << query << "\n" << reference << "\n";
  //
  //  auto ref_length{static_cast<int>(reference.length())};
  //
  //  StripedSmithWaterman::Aligner m_aligner{StripedSmithWaterman::Aligner{2, 5, 8, 6}};
  //  StripedSmithWaterman::Filter m_filter{StripedSmithWaterman::Filter{}};
  //  StripedSmithWaterman::Alignment m_alignment{};
  //
  //  m_aligner.Align(query.c_str(), reference.c_str(), ref_length, m_filter, &m_alignment,
  //                  ref_length / 2);
  //
  //  StripedSmithWaterman::print_alignment(query, reference, m_alignment);
  //
  std::cout << "sr: " << sr << "\n";

  return 0;
}
