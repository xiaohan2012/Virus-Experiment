from pymol import cmd
import os


#--------------------------------------VARS------------------------------------
InteractionName = ['any_i', 'backb', 'sidec', 'polar', 'hphob', 'hbacc', 'hbdon', 'aroma', 'charg']
InteractionButton = ['F1','F2','F3','F4','F5','F6','F7','F8','F9']
button_list = []
curSel=[0,0,0]

any_i=[1.00,1.00,1.00]
backb=[1.00,0.00,0.00]
sidec=[0.90,1.00,0.00]
polar=[0.00,0.00,1.00]
hphob=[0.25,0.25,0.00]
hbacc=[0.00,0.25,0.25]
hbdon=[0.25,0.00,0.25]
aroma=[0.20,0.40,0.50]
charg=[0.50,0.10,0.10]


#==============================================================================

def read_siftp(filename):
  
  try:
    #first line in file contains averaged sift

    with open(filename, 'r') as f:
        siftline =  f.readline()
    
    min_aa, sift = siftline.split(':')
    sift = sift.split(' ')
    print sift
    
    for interaction in range(9):
      selection = []
      for index in range(len(sift[interaction::9])):
          if index == len(sift[interaction::9]) - 1:
              print len(sift[interaction::9])
              print 'index',index,
              print 'interaction',sift[interaction::9][index]
          if sift[interaction::9][index].strip() != '0':
              if interaction == 8:
                  print sift[interaction::9][index]
                  print int(min_aa) + index
                  print sift[index:index + 9]
              selection.append(str(int(min_aa)+index))
      if len(selection) > 0:
          cmd.select ( InteractionButton[interaction]+'_'+InteractionName[interaction], 'resi ' + "+".join(selection) )
      if interaction == 0:
          cmd.show ( 'sticks', InteractionButton[interaction]+'_'+InteractionName[interaction] )
    
      cmd.remove ( 'hydro' )
      print InteractionButton[interaction]+'_'+InteractionName[interaction], 'resi ' + "+".join(selection) 
      del selection

    with open(filename) as f:
        f.readline()
        #print 'here we go'
        for line in f.readlines():
            try:
                #print line
                resi= line.split()[0]
                fp = line.split()[1:]
                fp =  ''.join(fp)
                #print resi,fp
                if fp != '000000000':
                    #print 'yes',fp
                    cmd.select('residue '+resi,'resi '+ resi)
                    cmd.show('sticks','resi '+ resi)
                    cmd.label('residue '+resi,'resi')
                else:
                    #cmd.select('residue '+resi,'resi '+ resi)
                    #cmd.label('residue '+resi,'"no"')
                    #print 'no',fp
                    pass
            except:
                pass
    
  except Exception, e:
    print "Failed to open file %s\n%s" % (filename, e)

#==============================================================================


def add_color(col2, butt):
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
  cmd.color('siff', InteractionButton[InteractionName.index(butt)]+'_'+butt )
  butt=[]

cmd.extend('read_siftp', read_siftp)

cmd.set_key( 'F1' , add_color, (any_i, 'any_i' ) )
cmd.set_key( 'F2' , add_color, (backb, 'backb' ) )
cmd.set_key( 'F3' , add_color, (sidec, 'sidec' ) )
cmd.set_key( 'F4' , add_color, (polar, 'polar' ) )
cmd.set_key( 'F5' , add_color, (hphob, 'hphob' ) )
cmd.set_key( 'F6' , add_color, (hbacc, 'hbacc' ) )
cmd.set_key( 'F7' , add_color, (hbdon, 'hbdon' ) )
cmd.set_key( 'F8' , add_color, (aroma, 'aroma' ) )
cmd.set_key( 'F9' , add_color, (charg, 'charg' ) )

