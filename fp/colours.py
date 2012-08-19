from pymol import cmd
import os


#--------------------------------------VARS------------------------------------
InteractionName = ['any_i', 'backb', 'sidec', 'polar', 'hphob', 'hbacc', 'hbdon', 'aroma', 'charg','env_polar','env_hphob','env_aroma',  'env_charg']
InteractionButton = ['f1','f2','f3','f4','f5','f6','f7','f8','f9','f10','f11','f12','f12']
button_list = []
curSel=[0,0,0]

any_i_exists=[1.00,1.00,1.00]
backb_exists=[1.00,0.00,0.00]
sidec_exists=[0.90,1.00,0.00]
polar_exists=[0.00,0.00,1.00]
hphob_exists=[0.25,0.25,0.00]
hbacc_exists=[0.00,0.25,0.25]
hbdon_exists=[0.25,0.00,0.25]
aroma_exists=[0.20,0.40,0.50]
charg_exists=[0.50,0.10,0.10]
env_polar_exists=[0.10,0.50,0.10]
env_hydrophobic_exists=[0.10,0.10,0.50]
env_aromatic_exists=[0.20,0.50,0.40]
env_charged_exists=[0.40,0.20,0.50]

any_i_strong=[1.00,1.00,1.00]
backb_strong=[1.00,0.00,0.00]
sidec_strong=[0.90,1.00,0.00]
polar_strong=[0.00,0.00,1.00]
hphob_strong=[0.25,0.25,0.00]
hbacc_strong=[0.00,0.25,0.25]
hbdon_strong=[0.25,0.00,0.25]
aroma_strong=[0.20,0.40,0.50]
charg_strong=[0.50,0.10,0.10]
env_polar_strong=[0.10,0.50,0.10]
env_hydrophobic_strong=[0.10,0.10,0.50]
env_aromatic_strong=[0.20,0.50,0.40]
env_charged_strong=[0.40,0.20,0.50]

#==============================================================================

def read_siftp(filename):
  
  try:
    #first line in file contains averaged sift

    with open(filename, 'r') as f:
        siftline =  f.readline()
    
    min_aa, sift = siftline.split(':')
    sift = sift.split(' ')
    print sift
    
    for interaction in range(13):
      selection_strong = []
      selection_exists = []
      for index in range(len(sift[interaction::13])):
          if sift[interaction::13][index].strip() == '10':
              selection_exists.append(str(int(min_aa)+index))
          elif sift[interaction::13][index].strip() == '11':
              selection_strong.append(str(int(min_aa)+index))

      if interaction == 0:#no distinction between strong and exists in ANY CONTACT type
          cmd.select ( InteractionName[interaction], 'resi ' + "+".join(selection_strong + selection_exists) )
          continue
          #cmd.show ( 'sticks', InteractionButton[interaction]+'_'+InteractionName[interaction] )
          
      if len(selection_strong) > 0:
          #cmd.select ( InteractionButton[interaction]+'_strong_'+InteractionName[interaction], 'resi ' + "+".join(selection_strong) )
          cmd.select ( InteractionName[interaction] + '_strong', 'resi ' + "+".join(selection_strong) )

      if len(selection_exists) > 0:
          #cmd.select ( InteractionButton[interaction]+'_exists_'+InteractionName[interaction], 'resi ' + "+".join(selection_exists) )
          cmd.select ( InteractionName[interaction] + '_exists', 'resi ' + "+".join(selection_exists) )
    
      cmd.remove ( 'hydro' )
      del selection_strong
      del selection_exists

  except Exception, e:
    print "Failed to open file %s\n%s" % (filename, e)

#==============================================================================


def add_color(col2, butt, t_):
# checking if the button was already hit
# if so, toggling the selection off

  button_list.append(butt)

  if button_list.count(butt) > 1:

    button_list.remove(butt)
    button_list.remove(butt)	
    cmd.util.cbag ('all')

  elif butt != button_list[0]:
    button_list.remove(button_list[0])
    cmd.util.cbag ('all')
    curSel = col2

  else:
    cmd.util.cbag ('all')
    curSel = col2

# colouring selection with specified colour scheme
  cmd.set_color('siff', curSel )
  #cmd.color('siff', InteractionButton[InteractionName.index(butt)]+'_' + t_ + '_'+butt )
  cmd.color('siff', butt )
  butt=[]

cmd.extend('read_siftp', read_siftp)
cmd.set_key( 'F1' , add_color, (any_i_exists, 'any_i_exists' ) )
cmd.set_key( 'F2' , add_color, (backb_exists, 'backb_exists' ) )
cmd.set_key( 'F3' , add_color, (sidec_exists, 'sidec_exists' ) )
cmd.set_key( 'F4' , add_color, (polar_exists, 'polar_exists' ) )
cmd.set_key( 'F5' , add_color, (hphob_exists, 'hphob_exists' ) )
cmd.set_key( 'F6' , add_color, (hbacc_exists, 'hbacc_exists' ) )
cmd.set_key( 'F7' , add_color, (hbdon_exists, 'hbdon_exists' ) )
cmd.set_key( 'F8' , add_color, (aroma_exists, 'aroma_exists' ) )
cmd.set_key( 'F9' , add_color, (charg_exists, 'charg_exists' ) )
cmd.set_key( 'F10' , add_color, (charg_exists, 'env_polar_exists' ) )
cmd.set_key( 'F11' , add_color, (charg_exists, 'env_hydrophobic_exists' ) )
cmd.set_key( 'F12' , add_color, (charg_exists, 'env_aromatic_exists' ) )
cmd.set_key( 'left' , add_color, (charg_exists, 'env_charged_exists' ) )

cmd.set_key( 'F1' , add_color, (any_i_strong, 'any_i_strong' ) )
cmd.set_key( 'F2' , add_color, (backb_strong, 'backb_strong' ) )
cmd.set_key( 'F3' , add_color, (sidec_strong, 'sidec_strong' ) )
cmd.set_key( 'F4' , add_color, (polar_strong, 'polar_strong' ) )
cmd.set_key( 'F5' , add_color, (hphob_strong, 'hphob_strong' ) )
cmd.set_key( 'F6' , add_color, (hbacc_strong, 'hbacc_strong' ) )
cmd.set_key( 'F7' , add_color, (hbdon_strong, 'hbdon_strong' ) )
cmd.set_key( 'F8' , add_color, (aroma_strong, 'aroma_strong' ) )
cmd.set_key( 'F9' , add_color, (charg_strong, 'charg_strong' ) )
cmd.set_key( 'F10' , add_color, (charg_strong, 'env_polar_strong' ) )
cmd.set_key( 'F11' , add_color, (charg_strong, 'env_hydrophobic_strong' ) )
cmd.set_key( 'F12' , add_color, (charg_strong, 'env_aromatic_strong' ) )
cmd.set_key( 'left' , add_color, (charg_strong, 'env_charged_strong' ) )
