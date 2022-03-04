#!/usr/bin/env python

import json
import glob


# def get_modules_of_genome_under_study(genome):
def count_nested_lists(l):
   count = 0 
   for e in l: 
      if isinstance(e, list):
         count = count + 1 + count_nested_lists(e)
   return count

def extract(lists, d):

   output = []

   if d == 1:
      return output.extend(lists)

   for sub_list in lists:
      extract(sub_list, d - 1)

def check_for_minus(alternative):

   complex_check = False

   if isinstance(alternative, list):
      tmp = []
      for i in alternative:    

         if "-" in i:
            tmp.append(i.split("-")[0])

         # elif "+" in i:
         #    split         = i.split("+")
         #    complex_check = True

         #    for j in split:
         #       tmp.append(j)
         else:
            tmp.append(i)
      alternative = tmp

   else: 
      if "-" in alternative: 
         alternative = alternative.split("-")[0]

      # if "+" in alternative: 
      #    alternative = alternative.split("+")

   
   return alternative



def main():
  
   f      = open("../../ref-dbs/module_definition_map.json", "r")
   mo_map = json.load(f)


   ncbi_id = str(435)
   for file in glob.glob("../../ref-dbs/kegg_genomes/" + ncbi_id + "/*.json"):
      with open(file, 'r') as f:
         genome_mos = json.load(f)

   # Species in-need of by-products
   count_of_species_modules = 0
   missing_steps_from_module_under_study = {}


   # Parse genome's KEGG modules and terms as found in the KEGG ORGANISMS db 
   for module, kos_on_its_own in genome_mos.items(): 

      list_of_kos_present                           = [k[3:] for k in kos_on_its_own.values()]
      count_of_species_modules                     += 1
      missing_steps_from_module_under_study[module] = {}

      if module != "md:M00376": 
         continue

      print("ON MY OWN: ", kos_on_its_own, "\n")
      print("The complete module: ", mo_map[module]['steps'], "\n")

      for index, step in enumerate(mo_map[module]['steps']):

         step_complet = False
         nested = count_nested_lists(step)

         """
         We will have a list (cases) for the various alternatives that will
         come up as sets of KOs that if found can complete the step
         We will only add cases in the missing_steps_from_module_under_study 
         dictionary if and oly if there is not already a way for a complete step
         """

         missing_steps_from_module_under_study[module][index] = []
         cases = []


         if nested > 0: 
            """
            In this case the children of the step, are the alternatives 
            and for each alternative, all the KOs are needed for the step
            to be complete
            """

            print("STEP: ", step)
            print("KOs present: ", list_of_kos_present)

            step_alternatives_map = {}

            for option, alternative in enumerate(step): 


               print("Alternative: >>", alternative)
               if isinstance(alternative, str):
                  alternative = [alternative]

               # Check for minus or/and pluses in step
               alternative_to_check = check_for_minus(alternative)
               checking = alternative_to_check[:]

               print("Alternative_to_check: ", alternative_to_check)

               # Investigate the alternatives
               for ko in alternative_to_check: 
                  print("IN THE LOOP: ", alternative_to_check)
                  print(ko)   
                  if ko in list_of_kos_present:
                     checking.remove(ko)

               if len(checking) == 0: 

                  print("\n\n >>>> ALL THE STEP IS HERE !!\n\n")

                  step_complet = True


               else: 
                  step_alternatives_map[option] = checking

            if step_complet == False: 
               
               missing_steps_from_module_under_study[module][index] = step_alternatives_map


         else:
            """
            In that case the alternatives, if any, are at the top level 
            and only on of them is needed for the step to be complete
            """

            step_without_signs = check_for_minus(step)
            step = step_without_signs[:]
            print("This is step: ", step)

            if isinstance(step, str):
               """
               The step can be as a string
               """
               if "+" in step: 

                  complex       = step.split("+")
                  complex_check = complex[:]
                  print(complex) is missi

                  for ko in complex:
                     if ko in list_of_kos_present:
                        complex_check.remove(ko)

                  if len(complex_check) == 0:
                     print("Complex complete on the species")

                  else:
                     missing_steps_from_module_under_study[module][index] = complex_check


            else:
               """
               Or as list
               """

               print("STEp:", step)

               for alternative in step:
                  """

                  """

                  if "+" in alternative: 

                     """
                     This is the case where a complex is on the step
                     """

                     print("Here is a complex:")

                     complex       = alternative.split("+")
                     complex_check = complex[:]
                     print(complex)

                     for ko in complex:
                        if ko in list_of_kos_present:
                           complex_check.remove(ko)

                     if len(complex_check) == 0:
                        print("Complex complete on the species")

                     else:
                        print(">>>> An incomplete complex!")
                        print(complex_check) 
                        missing_steps_from_module_under_study[module][index].append(complex_check)


                  else:

                     if alternative in list_of_kos_present:
                        break 

                     else:
                        missing_steps_from_module_under_study[module][index].append(alternative)

            print("\n\n~~~~\n")
            


            # if step_complet == False: 
            #    missing_steps_from_module_under_study[module][index] = missing


   missing_steps_from_module_under_study = {k: v for k, v in missing_steps_from_module_under_study.items() if v}
   print(" IN THE END:  ", missing_steps_from_module_under_study)

   return missing_steps_from_module_under_study





if __name__=="__main__":
   
   to_find = main()
   print(to_find)







# print("\n~~~~~\n")
# print(genome_mos['md:M00034'])
# print(missing_steps_from_module_under_study['md:M00034'])
# print("\n~~~~~\n")
# print("From its own: ", genome_mos['md:M00854'])
# print("To get from elsewher: ", missing_steps_from_module_under_study['md:M00854'])

