


def CreateAndOpen (Shot, Run, User, Device, Version, Subset, \
                Point, PrepareTrackPoint, PrepareTrackSubset):


    #!/usr/bin/env python

    #TODO:1    create 1 vtk (time 0, no scalars) and open in PV                         #DONE
    #TODO:2    create 1 vtk with one grid subset (time 0, no scalars), and open in PV   #DONE
    #TODO:3    create 2 blocks (.vtk) and open them in ParaView (geometry only, time 0) #DONE
    #TODO:4    create all the vtks and check all of them in ParaView                   #DONE
    #TODO:5    merge in a single vtkMultiBlock and open in PV (geo only and time 0)     #DONE
    #TODO:6    add electrons's scalars to those blocks and plot them in PV            #DONE
    #TODO:7    add all the scalars.                                                    #DONE
    #TODO:8    name all the scalars                                                    #DONE
    #TODO:9    deal with time. Good luck my friend                                      #DONE

    #TODO:10   create a list with indexes from "linear" datasets. Select them by hand   #DONE
    #TODO:11   create a list with time values (len)                                     #DONE 
    #TODO:12   create a list with the scalar's names corresponding to the indexes in 10 #DONE
    #TODO:13   create a list with the values of the scalars corresponding to 10         #DONE
    #TODO:14   write everything in a txt                                                #DONE
    #TODO:15   plot it in PV                                                            #DONE
    #TODO:16   add the name of the subset to the plot                                   #DONE
    #TODO:17   create main source                                                       #DONE
    #TODO:18   clean the code                                                           #DONE
    #TODO:19   plot the distance betwen the points and not the point IDs               #DONE





    from numpy import isnan

    import sys, os, shutil
    import imas
    import getopt
    import vtk

    class IDS_utilities:
        def __init__(self): 
            pass

        def ids_open(self, shot, run, user, device, version):
            #Open IMAS database (IDS)
            print("opening the IDS...")
            try:
                self.imas_obj = imas.ids(shot, run)
                self.imas_obj.open_env(user, device, version)
            except:
                print('Failed to open IDS. IDS does not exist!')
                self.state = False

            #we assume that our data is stored in the "edge_profiles" ggd
            self.imas_obj.edge_profiles.get()
            
            if self.imas_obj.isConnected():
                print('done!')
            else:
                print('Creation of data entry FAILED!')
                sys.exit()

        def ids_read_All(self, m=0, n=0, t=0):
            #Read IDS, get and set the necessary grid data (points geometry,
            #   connectivity array for 1D and 2D objects)
            
            print('reading the IDS...')

            points_geo = []
            list0D = []     # vertices
            list1D = []     # edges         
            list2D = []     # 2D cells /faces
            list3D = []     # 3D volumes

            
            try:
                #ggd = self.imas_obj.edge_profiles.ggd[t].grid
                ggd = self.imas_obj.edge_profiles.grid_ggd[t] #t = 0
                print("ggd opened!")
            except:
                print("Nothing found in edge profiles!"); sys.exit()


            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            #We will now get the number of objects for each dimension and the 
            #objects themselves

            #the number of different dimensions in our problem
            DIM = len(ggd.space[0].objects_per_dimension)


            if DIM >= 1:
                # 'Shortcut' variable to objects_per_dimension[0]
                ids_dim_0D = ggd.space[0].objects_per_dimension[0]
                num_obj_0D_all = len(ids_dim_0D.object)
            else:
                num_obj_0D_all = 0
            if DIM >= 2:
                ids_dim_1D = ggd.space[0].objects_per_dimension[1]
                num_obj_1D_all = len(ids_dim_1D.object)       
            else:
                num_obj_1D_all = 0
            if DIM >= 3:
                ids_dim_2D = ggd.space[0].objects_per_dimension[2]
                num_obj_2D_all = len(ids_dim_2D.object)
            else:
                num_obj_2D_all = 0
            if DIM >= 4:
                ids_dim_3D = ggd.space[0].objects_per_dimension[3]
                num_obj_3D_all = len(ids_dim_3D.object)
            else:
                num_obj_3D_all = 0
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            print("numer of points, ", num_obj_0D_all)
            print("numer of edges, ", num_obj_1D_all)
            print("numer of cells, ", num_obj_2D_all)


            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            #we get the points: their coordinates and their index 
            for i in range(num_obj_0D_all):
                points_geo.append(ids_dim_0D.object[i].geometry)
                list0D.append(int(ids_dim_0D.object[i].nodes[0]-1))
                #we add '-1' because IMAS is written in Fortran notation.

            #we get the indexes of the edges
            if num_obj_1D_all > 1:
                object = ids_dim_1D.object
                for i in range(num_obj_1D_all):
                    el1D = [int(object[i].nodes[j]-1) for j in range(2)]
                    #each line is formed by 2 nodes
                    list1D.append(el1D)
            #we get the indexes of the 2D cells
            if num_obj_2D_all > 1:
                object = ids_dim_2D.object
                m = len(object[0].nodes) #triangles, tetras... they have different 
                                        #number of nodes.
                for i in range(num_obj_2D_all):
                    el2D = [int(object[i].nodes[j]-1) for j in range(m)]
                    list2D.append(el2D)
                
            #we get the indexes of the 3D cells. In the most cases this line will
            #be useless.
            if num_obj_3D_all > 0:
                object = ids_dim_3D.object
                m = len(object[0].nodes)
                for i in range(num_obj_3D_all):
                    el3D = [int(object[i].nodes[j]-1) for j in range(m)]
                    list3D.append(el3D)
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


            return points_geo, list0D, list1D, list2D, list3D



        def ids_read_Scalars(self, t=0):

            ggd = self.imas_obj.edge_profiles.ggd[t]

            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            # these following variables will point directly to the IMAS "directory"
            # coontaining temperature, density, pressure...

            electrons = ggd.electrons           #there is only one type od electrons
            
            ions = []
            neutral = []

            for i in range(len(ggd.ion)):       #there can be multiple types of ions
                ions.append(ggd.ion[i])
            for i in range(len(neutral)):       #and the same for neutral
                neutral.append(ggd.neutral[i])
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            #these following variables contain e.g. 
            # [('temperature', <imas.ggd.etc.temp...>), (...)]
            # we need to define them in this way because there is no 
            # "temperature.name" or something like that. The first object is the 
            # name of the scalar and the second is the "path" to the values of the 
            # scalar.

            electrons_info =  [key for key in electrons.__dict__.items()]
            ions_info = []
            neutral_info = []

            for elem in ions:
                ions_info.append([key for key in elem.__dict__.items()])
            for elem in neutral:
                neutral_info.append([key for key in elem.__dict__.items()])
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
            # Now we will create one array for names and other for the "path"
            # to the scalars.
            names = []  # the names of the scalars (e.g. electrons_temperature) for
                        # ALL the fields (even the empty ones)

            paths = []  # the "paths" of the scalars for ALL the fields (temp, dens,
                        # ...) even the empty ones
                        
            #for more information, just print "names", "electrons_info", "paths"...
            names += ["electrons_"+elem[0] for elem in electrons_info]
            paths += [elem[1] for elem in electrons_info]

            for i in range(len(ions_info)):
                names += ["ions_"+str(i)+'_'+elem[0] for elem in ions_info[i]]
                paths += [elem[1] for elem in ions_info[i]]

            for i in range(len(neutral_info)):
                names += ["neutral_"+str(i)+'_'+elem[0] for elem in neutral_info[i]]
                paths += [elem[1] for elem in neutral_info[i]]
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
            #Now, we will "clean" both "names" and "paths", by creating the
            #following lists:
            
            SCALARS = []
            #    SCALARS -> [temperature[0].values, temperature[1].values, ...]
            #               [density[0].values, density[1].values, ...]
            #               [...]
            #               contains all scalars -> all values for each subset

            #e.g.
            #value = SCALARS[index corresponding to temp, dens, ...]
            #               [index corresponding to the subset (points, cells, ...)]
            #               [index of the scalar]

            #the temperature of the first point of the first subset (points) wil be
            #SCALARS[0][0][0]
            
            SC_INDEX = []
            #    SC_INDEX -> [subset_index related to temp[0], temp[1], temp[2]...]
            #                [subset_index related to dens[0], dens[1], dens[2]...]
            #                [...]

            NAMES = []


            #let's fill the lists
            for j, scalar in enumerate(paths):

                scalars = []    #[temperature[0].values, temperature[1].values, ...]
                                #... (nex iter)
        
                index = []      #[indexes of the matching grid subsets for temp.]
                                #...

                try:
                    for i in range(len(scalar)):  
                    #for each "path"


                        #note: some "0" are printed in some scalars, due to the way
                        #IMAS works. those "0" make the length of the "scalars" 
                        #array longer than the length of the cells to allocate those
                        #scallars, so we have to erase them.
                        #Also, we have found many 'nan' values in some scalars. We
                        #will also eliminate them.

                        #deal with random "0" and "nan" values
                        aux = [x for x in scalar[i].values if x != 0 and not isnan(x)]
                        if aux != []:
                            scalars.append(aux)
                            index.append(scalar[i].grid_subset_index)
                    

                except Exception as e:
                    #these errors are due to the fact that not every element of 
                    #"scalar" list is valid. To know the reason why we discard them,
                    #uncomment the following lines.

                    #print("some errors that we don't care about:")
                    #print(e)
                    pass

                SCALARS.append(scalars)
                SC_INDEX.append(index) 
                NAMES.append(names[j])
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            return SCALARS, SC_INDEX, NAMES


        def ids_close(self):
            #Close IDS database
            self.imas_obj.close()


        def ids_read_Subset(self, subset_index):
            #dimension 1 for points, 2 for edges and 3 for cells
            #reads the subset subset_index


            subset = \
                self.imas_obj.edge_profiles.grid_ggd[0].grid_subset[subset_index]
                #self.imas_obj.edge_profiles.ggd[0].grid.grid_subset[subset_index]
            
            name = subset.identifier.name #the name of the subset
            #subset_fortran_index = subset.identifier.index = subset_index + 1

            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            # those are the lists containing the cells that conform the subset. Each
            # subset will only have one type of cells. But just in case, let's keep
            # the three of them.
            index_0D = []
            index_1D = [] 
            index_2D = []

            for i in range(len(subset.element)):
                if subset.element[i].object[0].dimension == 1:
                    index_0D.append(subset.element[i].object[0].index-1)
                elif subset.element[i].object[0].dimension == 2:
                    index_1D.append(subset.element[i].object[0].index-1)
                elif subset.element[i].object[0].dimension == 3:
                    index_2D.append(subset.element[i].object[0].index-1)
            #we substract 1 because IMAS follows fortran notation.                
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      

            #we fix the same promlem with ceros that we had with the scalars.
            index_2D = [x for x in index_2D if x != -1]

            #print("there are "+str(len(index_0D)) +" points, "+str(len(index_1D)) \
            # + " lines   and "+str(len(index_2D)) +" cells in the g.s. \
            # "+str(subset_index))

            return index_0D, index_1D, index_2D, name

        def GetNumberOfSubsets(self):
            #return len(self.imas_obj.edge_profiles.ggd[0].grid.grid_subset)
            return len(self.imas_obj.edge_profiles.grid_ggd[0].grid_subset)
            

        def GetTimeArray(self):
            return [elem.time for elem in self.imas_obj.edge_profiles.ggd]

        def GetScIndexForNodes(self):
            for i in range(n):
                subset = \
                self.imas_obj.edge_profiles.grid_ggd[0].grid_subset[i]
                #self.imas_obj.edge_profiles.ggd[0].grid.grid_subset[i]

                name = subset.identifier.name #the name of the subset
                if name == "Nodes" or name == "nodes":
                    for j in range(len(SC_INDEX)):
                        if SC_INDEX[j] != []:
                            return SC_INDEX[j].index(i+1)




    def CreateVTK(index_0D, index_1D, index_2D, subset_index, SCALARS, \
                                            SC_INDEX, NAMES): 
        '''
        index_*D -> indices corresponding to verts, edges or lines
        subset_index -> the index of the subset, must match with scalar's one.
                        it must be in FORTRAN NOTATION.

        SCALARS ->  #[temperature[0].values, temperature[1].values, ...]
                    #[density[0].values, density[1].values, ...]
                    #[...]
                    contains all scalars -> all values for each subset
        e.g.
            value = SCALARS[index corresponding to temp, dens, ...]
                        [index corresponding to the subset (points, cells, ...)]
                        [index of the scalar]

            the temperature of the first point of the first subset (points) wil be
            SCALARS[0][0][0]

        SC_INDEX -> #[subset_index related to temp[0], temp[1], temp[2]...]
                    #[subset_index related to dens[0], dens[1], dens[2]...]
                    #[...]
        
        e.g.
            index = SC_INDEX[index corresponding to temp, dens, ...]
                            [index corresponding to the subset (points, cells, ...)]


        '''

        #the output of this function will be "datas", whitch is a object containing
        #all the grometry and the scalars selected in the input.
        datas = vtk.vtkPolyData()
    

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #First of all, we will add the grometry (points, cells...)

        #add points ---------------------------------------------  
        points = vtk.vtkPoints()
        if len(points_geo[0]) == 3:     #if 3D
            for elem in points_geo:
                points.InsertNextPoint(elem[0], elem[1], elem[2])   
        elif len(points_geo[0]) == 2:   #if 2D         
            for elem in points_geo:
                points.InsertNextPoint(elem[0], elem[1], 0)
        datas.SetPoints(points)
        
        #add 0D cells --------------------------------------------
        verts = vtk.vtkCellArray()
        for idx in index_0D:
            verts.InsertNextCell(1)
            verts.InsertCellPoint(idx)
        datas.SetVerts(verts)
        
        #add 1D cells --------------------------------------------
        lines = vtk.vtkCellArray()
        for i in index_1D:
            lines.InsertNextCell(2)
            lines.InsertCellPoint(cells1D[i][0])
            lines.InsertCellPoint(cells1D[i][1])
        datas.SetLines(lines)
        
        #add 2D cells --------------------------------------------
        if len(cells2D) > 0:
            tris = vtk.vtkCellArray()
            l = len(cells2D[0])

            for i in index_2D:
                tris.InsertNextCell(l)
                for j in range(l):
                    tris.InsertCellPoint(cells2D[i][j])

            datas.SetPolys(tris)  
        
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #Now, let's add all the scalars.


        VTKList = [] #this list will contain one vtkDoubleArray for each scalar, 
                    #which will be included later in "datas"
        
        


        for i in range(len(SCALARS)): 
        #for every scalar
            
            if subset_index in SC_INDEX[i]:
            #if the subset "subset_index" has some scalars associated

                if len(SCALARS[i]) > 0:   
                    #if those "associated scalars" actually exist
            
                    VTKList.append(vtk.vtkDoubleArray())
                    #every property will go to the vtkDoubleArray
                    
                    VTKList[-1].SetName(NAMES[i])

                    for elem in SCALARS[i][SC_INDEX[i].index(subset_index)]:   
                    #for each value (number) associated to subset_index and to the
                    #i th scalar
                        VTKList[-1].InsertNextValue(elem)

                else:
                    print("scalaes[i] EMPTY")

        #let's add the scalars to "datas"
        for elem in VTKList:
            datas.GetCellData().AddArray(elem)
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        return datas





    def sort(some_list):
        #['a10', 'a20', 'a5'] -> ['a5', 'a10', 'a20']
        numbers = []
        sortedd= [0] * len(some_list)

        for elem in some_list:
            aux = ''
            for char in elem:
                if char >='0' and char <='9':
                    aux += char
            if aux != '':
                numbers.append(int(aux))

        for i, elem in enumerate(numbers):
            sortedd[elem] = some_list[i]
        return sortedd


    shot = Shot     
    run = Run            
    user = User     
    device = Device          
    version = Version             
    point_ID = Point   
    ChosenSubset = Subset
    TrackSubset = PrepareTrackSubset
    TrackPoint = PrepareTrackPoint


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # if ChosenSubset == 999, a table with the available subsets will be displayed

    name_final = ['Inner target', 'Inner throat', 'Outer throat', 'Outer target', \
                    'Core cut', 'PFR cut', 'Inner PFR wall', 'Outer PFR wall', \
                    'Inner baffle', 'Outer baffle']

    if ChosenSubset == 999:
        space = ' '

        max_len = len(name_final[0])

        for elem in name_final[1:]:
            if len(elem) > max_len:
                max_len = len(elem)

        print("\n\n Available subsets:\n")      #print the formatted table
        print("|    NAME  "+(max_len-12)*' '+"   | SUBSET |")  
        print("|"+max_len*'-'+"+---------|")

        for i, elem in enumerate(name_final):

            sp_units = max_len-len(name_final[i])
            sp = space*sp_units
            print('|'+elem + sp + '     ' + str(i) + '    |')
        print("\n\n")

        return 0
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #we will clean the temp folder before writting new data
    folder = './temp_plugin' 
    if not 'temp_plugin' in os.listdir('.'):
        os.mkdir(folder)
    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # We will now inicialize the IDS_utilities class
    ids_utils = IDS_utilities()
    ids_utils.ids_open(shot, run, user, device, version)

    #...and get information about time
    timeArray = ids_utils.GetTimeArray()
    NumberOfTimeSlices = len(timeArray)

    #... and the information about the mesh from the first timeslice
    points_geo, obj_0D_list, cells1D, cells2D, cells3D \
                                    = ids_utils.ids_read_All(t = 0)

    #... and the number of subsets
    n = ids_utils.GetNumberOfSubsets()
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # this is the main bucle. We write here all the .vtm that are part of our 
    # animation

    for t in range(NumberOfTimeSlices):
        fileName = './temp_plugin/multiblock_'+str(t)+'.vtm'

        #for each t, we create a multiBlock
        multiBlock = vtk.vtkMultiBlockDataSet()
        multiBlock.SetNumberOfBlocks(n)
        SCALARS, SC_INDEX, NAMES = ids_utils.ids_read_Scalars(t = t)

        for i in range(n): #we write n .vtp for each multiblock 
            index_0D, index_1D, index_2D, name = ids_utils.ids_read_Subset(i)
            #TODO: read_Subset should be outside this bucle.

            #we write each vtp, n of them for each block. i+1 is the subset_index
            #written in fortran notation.
            datas = CreateVTK(index_0D, index_1D, index_2D, \
                                i+1, SCALARS, SC_INDEX, NAMES)

            multiBlock.SetBlock(i, datas)
            a = multiBlock.GetMetaData(i)
            a.Set(multiBlock.NAME(), name)

        #just a progress checker
        perc = (t+1)*100.0/NumberOfTimeSlices
        print("writing " + fileName + "; Total progress: %.1f  %%. \
        Current: timeStep %d / %d" %(perc, t+1, NumberOfTimeSlices))

        writer = vtk.vtkXMLMultiBlockDataWriter()
        writer.SetDataModeToAscii()

        writer.SetFileName(fileName)
        writer.SetInputData(multiBlock)
        writer.Write()
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




    print("done!\n\n")

    print("-----------------------------------------------------------------------")



    if TrackSubset:
        #if the user has decided to wtite the data for the SubsetTracker

        print("preparing the data for SubsetTracker!\n")

        #available subset_grids for PythonView ( lines ):
        #
        #   inner tatget    inner throat    outer throat
        #   outer target    core cut        PFR cut
        #   inner PRF wall  outer PRF wall  inner baffle  
        #   outer baffle
        #
        #available subset_grids for PythonView ( points ):
        #
        #   inner Midplane  Outer Midplane
        #
        #------------------------------------------------


        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #from now we work only with the selected subset_index. First we extract it 
        #from the whole list.

        #these lists will contain values for every subset. We only want one of them.
        index_1D_aux = []
        name_aux = []
        subset_index_aux = []

        for i in range(n):
            index_0D, index_1D, index_2D, name = ids_utils.ids_read_Subset(i)
            index_1D_aux.append(index_1D)
            name_aux.append(name)
            subset_index_aux.append(i+1)


        index_1D_final = index_1D_aux[name_aux.index(name_final[ChosenSubset])]
        #indexes to the lines (2 elem list) forming the subset
        #(same length as name_final)
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        OurScalars = [] #list containing all the values for all the scalars
            
            #OurScalars = [[temperature.values], [dens.values], ...] (time 0)
            #        and so on for further times
            #each OurScalars[i] contains len(POINTS) elements
            #e.g.: value = OurScalars[timestep][scalar_index]

        print("taking from the whole list only the values of all the scalars")
        print("related to the points of the choosen subset...\n")


        #let's fill OurScalars
        for t in range(NumberOfTimeSlices):
            scalar_names = []

            SCALARS, SC_INDEX, NAMES = ids_utils.ids_read_Scalars(t = t)
            aux = []    #we add the scalars

            for i, elem in enumerate(SCALARS):
                if len(elem) > 0:
                    scalar_names.append(NAMES[i])
                    aux.append(elem[ids_utils.GetScIndexForNodes()])
            OurScalars.append(aux)  #we add the times

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        chosen_name = name_final[ChosenSubset]
        #the name of the chosen subset

        print("you have chosen: " + chosen_name)


        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #we will add now the points of our subset.

        print("taking from the whole list the indexes of the points \
        forming the choosen subset...\n")

        #we will make now a list containing the indexes to the points that are part 
        #of the chosen subset
        indexToPoints = []  #we first add the first 2 nodes
        indexToPoints += cells1D[index_1D_final[0]]

        #now the rest
        for k in range(1, len(index_1D_final)):
            aux = cells1D[index_1D_final[k]][1]
            indexToPoints.append(aux)
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #we will add now the scalars

        scalar_values = []  
                        #the value of the scalars pointed by elem in indexToPoints.
                        #scalar_values = [[temperature.values related to the nodes],
                        #                [dens.values related to the nodes], 
                        #                ...]

        for t in range(NumberOfTimeSlices):
            aux = []
            for scalar in OurScalars[t]:
                aux.append([scalar[node] for node in indexToPoints])
            scalar_values.append(aux)
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #Now we're going to take the coordinates of the nodes and make a list from
        #their distance, starting in 0 and ending in the correspondent distance.
        #this will be a great improvement for the 3D plot.
        
        def distance2D(p1, p2):
        #we will calculate distances with this function.    
            return  ((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)**0.5

        distances = [0] #this will be our list (firs element in position 0)

        for i in range(1, len(indexToPoints)):
            distances.append(
                distance2D(points_geo[i], points_geo[i-1])+distances[-1])
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #we write everything now



        outputName = "./temp_plugin/tunnel_plot.txt"
        g = open(outputName,"w")


        #first some information about the mesh
        g.write(chosen_name)            #the name of the subset
        g.write('\n')
        for position in distances:     #list for "x" variable for the plot
            g.write(str(position) + ' ')
        g.write('\n')

        for time in timeArray:          #list for "time" varibale
            g.write(str(time) + ' ')
        g.write('\n')


        #now the titles: "time", "point_ID" and scalars
        g.write("time" + ' ')
        g.write("point_ID" + ' ')
        for name in scalar_names:
            g.write(name + ' ')
        g.write("\n")


        print("writing everything on "+outputName+"...\n")


        for t in range(NumberOfTimeSlices):             #for every timestep
            for i in range(len(indexToPoints)):         #for every point
                g.write(str(t)+' ')
                g.write(str(indexToPoints[i])+' ')
                for elem in scalar_values[t]:           #for every scalar
                    g.write(str(elem[i]))
                    g.write(" ")
                g.write('\n')

        g.close()

        print("done!\n\n")
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    print("-----------------------------------------------------------------------")



    if TrackPoint:  #if we have decided to track a single point
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #we get here all the needed scalars

        print("preparing the data for PointTracker...")
        outputName = "./temp_plugin/trackpoint.txt"
        f = open(outputName,"w")


        OurScalars = [] #value = OurScalars[timestep][scalar (temp, dens...)]
        #this list contains the values of the sccalars for each timestep
        for t in range(NumberOfTimeSlices):
            scalar_names = [] 

            SCALARS, SC_INDEX, NAMES = ids_utils.ids_read_Scalars(t = t) #first t = 0
            aux = []
            for i, elem in enumerate(SCALARS):
                if len(elem) > 0:
                    scalar_names.append(NAMES[i])
                    aux.append(elem[ids_utils.GetScIndexForNodes()][point_ID])
            OurScalars.append(aux)
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #we write everything now

        f.write("time" + ' ')
        for i in range(len(scalar_names)):
            f.write(scalar_names[i] + ' ')
        f.write("\n")

        print("writing everything on "+outputName+"...\n")

        for t in range(NumberOfTimeSlices):
            
            f.write(str(timeArray[t]))
            f.write(" ")

            for elem in OurScalars[t]:
                f.write(str(elem))
                f.write(" ")
            f.write('\n')   
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        f.close()


    ids_utils.ids_close()
    print("-----------------------------------------------------------------------")
    print("Finished!")


#uncomment if running on a ParaView's Programmable Source

shot = 22108  
run = 1000
user = 'g2mcarpi'                       
device = 'imas20'           
version = '3'              
PrepareTrackSubset = True      #if you want to write a file for tracking a subset
PrepareTrackPoint = True       #if you want to track a point   
Subset = 5              #the subset that you want to track
Point = 900          #the point you want to track

CreateAndOpen(shot, run, user, device, version, Subset, \
                Point, PrepareTrackPoint, PrepareTrackSubset)




def CreateAndOpen (Shot, Run, User, Device, Version, Subset, \
                Point, PrepareTrackPoint, PrepareTrackSubset):


    #!/usr/bin/env python

    #TODO:1    create 1 vtk (time 0, no scalars) and open in PV                         #DONE
    #TODO:2    create 1 vtk with one grid subset (time 0, no scalars), and open in PV   #DONE
    #TODO:3    create 2 blocks (.vtk) and open them in ParaView (geometry only, time 0) #DONE
    #TODO:4    create all the vtks and check all of them in ParaView                   #DONE
    #TODO:5    merge in a single vtkMultiBlock and open in PV (geo only and time 0)     #DONE
    #TODO:6    add electrons's scalars to those blocks and plot them in PV            #DONE
    #TODO:7    add all the scalars.                                                    #DONE
    #TODO:8    name all the scalars                                                    #DONE
    #TODO:9    deal with time. Good luck my friend                                      #DONE

    #TODO:10   create a list with indexes from "linear" datasets. Select them by hand   #DONE
    #TODO:11   create a list with time values (len)                                     #DONE 
    #TODO:12   create a list with the scalar's names corresponding to the indexes in 10 #DONE
    #TODO:13   create a list with the values of the scalars corresponding to 10         #DONE
    #TODO:14   write everything in a txt                                                #DONE
    #TODO:15   plot it in PV                                                            #DONE
    #TODO:16   add the name of the subset to the plot                                   #DONE
    #TODO:17   create main source                                                       #DONE
    #TODO:18   clean the code                                                           #DONE
    #TODO:19   plot the distance betwen the points and not the point IDs               #DONE





    from numpy import isnan

    import sys, os, shutil
    import imas
    import getopt
    import vtk

    class IDS_utilities:
        def __init__(self): 
            pass

        def ids_open(self, shot, run, user, device, version):
            #Open IMAS database (IDS)
            print("opening the IDS...")
            try:
                self.imas_obj = imas.ids(shot, run)
                self.imas_obj.open_env(user, device, version)
            except:
                print('Failed to open IDS. IDS does not exist!')
                self.state = False

            #we assume that our data is stored in the "edge_profiles" ggd
            self.imas_obj.edge_profiles.get()
            
            if self.imas_obj.isConnected():
                print('done!')
            else:
                print('Creation of data entry FAILED!')
                sys.exit()

        def ids_read_All(self, m=0, n=0, t=0):
            #Read IDS, get and set the necessary grid data (points geometry,
            #   connectivity array for 1D and 2D objects)
            
            print('reading the IDS...')

            points_geo = []
            list0D = []     # vertices
            list1D = []     # edges         
            list2D = []     # 2D cells /faces
            list3D = []     # 3D volumes

            
            try:
                #ggd = self.imas_obj.edge_profiles.ggd[t].grid
                ggd = self.imas_obj.edge_profiles.grid_ggd[t] #t = 0
                print("ggd opened!")
            except:
                print("Nothing found in edge profiles!"); sys.exit()


            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            #We will now get the number of objects for each dimension and the 
            #objects themselves

            #the number of different dimensions in our problem
            DIM = len(ggd.space[0].objects_per_dimension)


            if DIM >= 1:
                # 'Shortcut' variable to objects_per_dimension[0]
                ids_dim_0D = ggd.space[0].objects_per_dimension[0]
                num_obj_0D_all = len(ids_dim_0D.object)
            else:
                num_obj_0D_all = 0
            if DIM >= 2:
                ids_dim_1D = ggd.space[0].objects_per_dimension[1]
                num_obj_1D_all = len(ids_dim_1D.object)       
            else:
                num_obj_1D_all = 0
            if DIM >= 3:
                ids_dim_2D = ggd.space[0].objects_per_dimension[2]
                num_obj_2D_all = len(ids_dim_2D.object)
            else:
                num_obj_2D_all = 0
            if DIM >= 4:
                ids_dim_3D = ggd.space[0].objects_per_dimension[3]
                num_obj_3D_all = len(ids_dim_3D.object)
            else:
                num_obj_3D_all = 0
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            print("numer of points, ", num_obj_0D_all)
            print("numer of edges, ", num_obj_1D_all)
            print("numer of cells, ", num_obj_2D_all)


            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            #we get the points: their coordinates and their index 
            for i in range(num_obj_0D_all):
                points_geo.append(ids_dim_0D.object[i].geometry)
                list0D.append(int(ids_dim_0D.object[i].nodes[0]-1))
                #we add '-1' because IMAS is written in Fortran notation.

            #we get the indexes of the edges
            if num_obj_1D_all > 1:
                object = ids_dim_1D.object
                for i in range(num_obj_1D_all):
                    el1D = [int(object[i].nodes[j]-1) for j in range(2)]
                    #each line is formed by 2 nodes
                    list1D.append(el1D)
            #we get the indexes of the 2D cells
            if num_obj_2D_all > 1:
                object = ids_dim_2D.object
                m = len(object[0].nodes) #triangles, tetras... they have different 
                                        #number of nodes.
                for i in range(num_obj_2D_all):
                    el2D = [int(object[i].nodes[j]-1) for j in range(m)]
                    list2D.append(el2D)
                
            #we get the indexes of the 3D cells. In the most cases this line will
            #be useless.
            if num_obj_3D_all > 0:
                object = ids_dim_3D.object
                m = len(object[0].nodes)
                for i in range(num_obj_3D_all):
                    el3D = [int(object[i].nodes[j]-1) for j in range(m)]
                    list3D.append(el3D)
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


            return points_geo, list0D, list1D, list2D, list3D



        def ids_read_Scalars(self, t=0):

            ggd = self.imas_obj.edge_profiles.ggd[t]

            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            # these following variables will point directly to the IMAS "directory"
            # coontaining temperature, density, pressure...

            electrons = ggd.electrons           #there is only one type od electrons
            
            ions = []
            neutral = []

            for i in range(len(ggd.ion)):       #there can be multiple types of ions
                ions.append(ggd.ion[i])
            for i in range(len(neutral)):       #and the same for neutral
                neutral.append(ggd.neutral[i])
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            #these following variables contain e.g. 
            # [('temperature', <imas.ggd.etc.temp...>), (...)]
            # we need to define them in this way because there is no 
            # "temperature.name" or something like that. The first object is the 
            # name of the scalar and the second is the "path" to the values of the 
            # scalar.

            electrons_info =  [key for key in electrons.__dict__.items()]
            ions_info = []
            neutral_info = []

            for elem in ions:
                ions_info.append([key for key in elem.__dict__.items()])
            for elem in neutral:
                neutral_info.append([key for key in elem.__dict__.items()])
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
            # Now we will create one array for names and other for the "path"
            # to the scalars.
            names = []  # the names of the scalars (e.g. electrons_temperature) for
                        # ALL the fields (even the empty ones)

            paths = []  # the "paths" of the scalars for ALL the fields (temp, dens,
                        # ...) even the empty ones
                        
            #for more information, just print "names", "electrons_info", "paths"...
            names += ["electrons_"+elem[0] for elem in electrons_info]
            paths += [elem[1] for elem in electrons_info]

            for i in range(len(ions_info)):
                names += ["ions_"+str(i)+'_'+elem[0] for elem in ions_info[i]]
                paths += [elem[1] for elem in ions_info[i]]

            for i in range(len(neutral_info)):
                names += ["neutral_"+str(i)+'_'+elem[0] for elem in neutral_info[i]]
                paths += [elem[1] for elem in neutral_info[i]]
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
            #Now, we will "clean" both "names" and "paths", by creating the
            #following lists:
            
            SCALARS = []
            #    SCALARS -> [temperature[0].values, temperature[1].values, ...]
            #               [density[0].values, density[1].values, ...]
            #               [...]
            #               contains all scalars -> all values for each subset

            #e.g.
            #value = SCALARS[index corresponding to temp, dens, ...]
            #               [index corresponding to the subset (points, cells, ...)]
            #               [index of the scalar]

            #the temperature of the first point of the first subset (points) wil be
            #SCALARS[0][0][0]
            
            SC_INDEX = []
            #    SC_INDEX -> [subset_index related to temp[0], temp[1], temp[2]...]
            #                [subset_index related to dens[0], dens[1], dens[2]...]
            #                [...]

            NAMES = []


            #let's fill the lists
            for j, scalar in enumerate(paths):

                scalars = []    #[temperature[0].values, temperature[1].values, ...]
                                #... (nex iter)
        
                index = []      #[indexes of the matching grid subsets for temp.]
                                #...

                try:
                    for i in range(len(scalar)):  
                    #for each "path"


                        #note: some "0" are printed in some scalars, due to the way
                        #IMAS works. those "0" make the length of the "scalars" 
                        #array longer than the length of the cells to allocate those
                        #scallars, so we have to erase them.
                        #Also, we have found many 'nan' values in some scalars. We
                        #will also eliminate them.

                        #deal with random "0" and "nan" values
                        aux = [x for x in scalar[i].values if x != 0 and not isnan(x)]
                        if aux != []:
                            scalars.append(aux)
                            index.append(scalar[i].grid_subset_index)
                    

                except Exception as e:
                    #these errors are due to the fact that not every element of 
                    #"scalar" list is valid. To know the reason why we discard them,
                    #uncomment the following lines.

                    #print("some errors that we don't care about:")
                    #print(e)
                    pass

                SCALARS.append(scalars)
                SC_INDEX.append(index) 
                NAMES.append(names[j])
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            return SCALARS, SC_INDEX, NAMES


        def ids_close(self):
            #Close IDS database
            self.imas_obj.close()


        def ids_read_Subset(self, subset_index):
            #dimension 1 for points, 2 for edges and 3 for cells
            #reads the subset subset_index


            subset = \
                self.imas_obj.edge_profiles.grid_ggd[0].grid_subset[subset_index]
                #self.imas_obj.edge_profiles.ggd[0].grid.grid_subset[subset_index]
            
            name = subset.identifier.name #the name of the subset
            #subset_fortran_index = subset.identifier.index = subset_index + 1

            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            # those are the lists containing the cells that conform the subset. Each
            # subset will only have one type of cells. But just in case, let's keep
            # the three of them.
            index_0D = []
            index_1D = [] 
            index_2D = []

            for i in range(len(subset.element)):
                if subset.element[i].object[0].dimension == 1:
                    index_0D.append(subset.element[i].object[0].index-1)
                elif subset.element[i].object[0].dimension == 2:
                    index_1D.append(subset.element[i].object[0].index-1)
                elif subset.element[i].object[0].dimension == 3:
                    index_2D.append(subset.element[i].object[0].index-1)
            #we substract 1 because IMAS follows fortran notation.                
            #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      

            #we fix the same promlem with ceros that we had with the scalars.
            index_2D = [x for x in index_2D if x != -1]

            #print("there are "+str(len(index_0D)) +" points, "+str(len(index_1D)) \
            # + " lines   and "+str(len(index_2D)) +" cells in the g.s. \
            # "+str(subset_index))

            return index_0D, index_1D, index_2D, name

        def GetNumberOfSubsets(self):
            #return len(self.imas_obj.edge_profiles.ggd[0].grid.grid_subset)
            return len(self.imas_obj.edge_profiles.grid_ggd[0].grid_subset)
            

        def GetTimeArray(self):
            return [elem.time for elem in self.imas_obj.edge_profiles.ggd]

        def GetScIndexForNodes(self):
            for i in range(n):
                subset = \
                self.imas_obj.edge_profiles.grid_ggd[0].grid_subset[i]
                #self.imas_obj.edge_profiles.ggd[0].grid.grid_subset[i]

                name = subset.identifier.name #the name of the subset
                if name == "Nodes" or name == "nodes":
                    for j in range(len(SC_INDEX)):
                        if SC_INDEX[j] != []:
                            return SC_INDEX[j].index(i+1)




    def CreateVTK(index_0D, index_1D, index_2D, subset_index, SCALARS, \
                                            SC_INDEX, NAMES): 
        '''
        index_*D -> indices corresponding to verts, edges or lines
        subset_index -> the index of the subset, must match with scalar's one.
                        it must be in FORTRAN NOTATION.

        SCALARS ->  #[temperature[0].values, temperature[1].values, ...]
                    #[density[0].values, density[1].values, ...]
                    #[...]
                    contains all scalars -> all values for each subset
        e.g.
            value = SCALARS[index corresponding to temp, dens, ...]
                        [index corresponding to the subset (points, cells, ...)]
                        [index of the scalar]

            the temperature of the first point of the first subset (points) wil be
            SCALARS[0][0][0]

        SC_INDEX -> #[subset_index related to temp[0], temp[1], temp[2]...]
                    #[subset_index related to dens[0], dens[1], dens[2]...]
                    #[...]
        
        e.g.
            index = SC_INDEX[index corresponding to temp, dens, ...]
                            [index corresponding to the subset (points, cells, ...)]


        '''

        #the output of this function will be "datas", whitch is a object containing
        #all the grometry and the scalars selected in the input.
        datas = vtk.vtkPolyData()
    

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #First of all, we will add the grometry (points, cells...)

        #add points ---------------------------------------------  
        points = vtk.vtkPoints()
        if len(points_geo[0]) == 3:     #if 3D
            for elem in points_geo:
                points.InsertNextPoint(elem[0], elem[1], elem[2])   
        elif len(points_geo[0]) == 2:   #if 2D         
            for elem in points_geo:
                points.InsertNextPoint(elem[0], elem[1], 0)
        datas.SetPoints(points)
        
        #add 0D cells --------------------------------------------
        verts = vtk.vtkCellArray()
        for idx in index_0D:
            verts.InsertNextCell(1)
            verts.InsertCellPoint(idx)
        datas.SetVerts(verts)
        
        #add 1D cells --------------------------------------------
        lines = vtk.vtkCellArray()
        for i in index_1D:
            lines.InsertNextCell(2)
            lines.InsertCellPoint(cells1D[i][0])
            lines.InsertCellPoint(cells1D[i][1])
        datas.SetLines(lines)
        
        #add 2D cells --------------------------------------------
        if len(cells2D) > 0:
            tris = vtk.vtkCellArray()
            l = len(cells2D[0])

            for i in index_2D:
                tris.InsertNextCell(l)
                for j in range(l):
                    tris.InsertCellPoint(cells2D[i][j])

            datas.SetPolys(tris)  
        
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #Now, let's add all the scalars.


        VTKList = [] #this list will contain one vtkDoubleArray for each scalar, 
                    #which will be included later in "datas"
        
        


        for i in range(len(SCALARS)): 
        #for every scalar
            
            if subset_index in SC_INDEX[i]:
            #if the subset "subset_index" has some scalars associated

                if len(SCALARS[i]) > 0:   
                    #if those "associated scalars" actually exist
            
                    VTKList.append(vtk.vtkDoubleArray())
                    #every property will go to the vtkDoubleArray
                    
                    VTKList[-1].SetName(NAMES[i])

                    for elem in SCALARS[i][SC_INDEX[i].index(subset_index)]:   
                    #for each value (number) associated to subset_index and to the
                    #i th scalar
                        VTKList[-1].InsertNextValue(elem)

                else:
                    print("scalaes[i] EMPTY")

        #let's add the scalars to "datas"
        for elem in VTKList:
            datas.GetCellData().AddArray(elem)
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        return datas





    def sort(some_list):
        #['a10', 'a20', 'a5'] -> ['a5', 'a10', 'a20']
        numbers = []
        sortedd= [0] * len(some_list)

        for elem in some_list:
            aux = ''
            for char in elem:
                if char >='0' and char <='9':
                    aux += char
            if aux != '':
                numbers.append(int(aux))

        for i, elem in enumerate(numbers):
            sortedd[elem] = some_list[i]
        return sortedd


    shot = Shot     
    run = Run            
    user = User     
    device = Device          
    version = Version             
    point_ID = Point   
    ChosenSubset = Subset
    TrackSubset = PrepareTrackSubset
    TrackPoint = PrepareTrackPoint


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # if ChosenSubset == 999, a table with the available subsets will be displayed

    name_final = ['Inner target', 'Inner throat', 'Outer throat', 'Outer target', \
                    'Core cut', 'PFR cut', 'Inner PFR wall', 'Outer PFR wall', \
                    'Inner baffle', 'Outer baffle']

    if ChosenSubset == 999:
        space = ' '

        max_len = len(name_final[0])

        for elem in name_final[1:]:
            if len(elem) > max_len:
                max_len = len(elem)

        print("\n\n Available subsets:\n")      #print the formatted table
        print("|    NAME  "+(max_len-12)*' '+"   | SUBSET |")  
        print("|"+max_len*'-'+"+---------|")

        for i, elem in enumerate(name_final):

            sp_units = max_len-len(name_final[i])
            sp = space*sp_units
            print('|'+elem + sp + '     ' + str(i) + '    |')
        print("\n\n")

        return 0
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #we will clean the temp folder before writting new data
    folder = './temp_plugin' 
    if not 'temp_plugin' in os.listdir('.'):
        os.mkdir(folder)
    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # We will now inicialize the IDS_utilities class
    ids_utils = IDS_utilities()
    ids_utils.ids_open(shot, run, user, device, version)

    #...and get information about time
    timeArray = ids_utils.GetTimeArray()
    NumberOfTimeSlices = len(timeArray)

    #... and the information about the mesh from the first timeslice
    points_geo, obj_0D_list, cells1D, cells2D, cells3D \
                                    = ids_utils.ids_read_All(t = 0)

    #... and the number of subsets
    n = ids_utils.GetNumberOfSubsets()
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # this is the main bucle. We write here all the .vtm that are part of our 
    # animation

    for t in range(NumberOfTimeSlices):
        fileName = './temp_plugin/multiblock_'+str(t)+'.vtm'

        #for each t, we create a multiBlock
        multiBlock = vtk.vtkMultiBlockDataSet()
        multiBlock.SetNumberOfBlocks(n)
        SCALARS, SC_INDEX, NAMES = ids_utils.ids_read_Scalars(t = t)

        for i in range(n): #we write n .vtp for each multiblock 
            index_0D, index_1D, index_2D, name = ids_utils.ids_read_Subset(i)
            #TODO: read_Subset should be outside this bucle.

            #we write each vtp, n of them for each block. i+1 is the subset_index
            #written in fortran notation.
            datas = CreateVTK(index_0D, index_1D, index_2D, \
                                i+1, SCALARS, SC_INDEX, NAMES)

            multiBlock.SetBlock(i, datas)
            a = multiBlock.GetMetaData(i)
            a.Set(multiBlock.NAME(), name)

        #just a progress checker
        perc = (t+1)*100.0/NumberOfTimeSlices
        print("writing " + fileName + "; Total progress: %.1f  %%. \
        Current: timeStep %d / %d" %(perc, t+1, NumberOfTimeSlices))

        writer = vtk.vtkXMLMultiBlockDataWriter()
        writer.SetDataModeToAscii()

        writer.SetFileName(fileName)
        writer.SetInputData(multiBlock)
        writer.Write()
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




    print("done!\n\n")

    print("-----------------------------------------------------------------------")



    if TrackSubset:
        #if the user has decided to wtite the data for the SubsetTracker

        print("preparing the data for SubsetTracker!\n")

        #available subset_grids for PythonView ( lines ):
        #
        #   inner tatget    inner throat    outer throat
        #   outer target    core cut        PFR cut
        #   inner PRF wall  outer PRF wall  inner baffle  
        #   outer baffle
        #
        #available subset_grids for PythonView ( points ):
        #
        #   inner Midplane  Outer Midplane
        #
        #------------------------------------------------


        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #from now we work only with the selected subset_index. First we extract it 
        #from the whole list.

        #these lists will contain values for every subset. We only want one of them.
        index_1D_aux = []
        name_aux = []
        subset_index_aux = []

        for i in range(n):
            index_0D, index_1D, index_2D, name = ids_utils.ids_read_Subset(i)
            index_1D_aux.append(index_1D)
            name_aux.append(name)
            subset_index_aux.append(i+1)


        index_1D_final = index_1D_aux[name_aux.index(name_final[ChosenSubset])]
        #indexes to the lines (2 elem list) forming the subset
        #(same length as name_final)
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        OurScalars = [] #list containing all the values for all the scalars
            
            #OurScalars = [[temperature.values], [dens.values], ...] (time 0)
            #        and so on for further times
            #each OurScalars[i] contains len(POINTS) elements
            #e.g.: value = OurScalars[timestep][scalar_index]

        print("taking from the whole list only the values of all the scalars")
        print("related to the points of the choosen subset...\n")


        #let's fill OurScalars
        for t in range(NumberOfTimeSlices):
            scalar_names = []

            SCALARS, SC_INDEX, NAMES = ids_utils.ids_read_Scalars(t = t)
            aux = []    #we add the scalars

            for i, elem in enumerate(SCALARS):
                if len(elem) > 0:
                    scalar_names.append(NAMES[i])
                    aux.append(elem[ids_utils.GetScIndexForNodes()])
            OurScalars.append(aux)  #we add the times

        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        chosen_name = name_final[ChosenSubset]
        #the name of the chosen subset

        print("you have chosen: " + chosen_name)


        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #we will add now the points of our subset.

        print("taking from the whole list the indexes of the points \
corresponding the choosen subset...\n")

        #we will make now a list containing the indexes to the points that are part 
        #of the chosen subset
        indexToPoints = []  #we first add the first 2 nodes
        indexToPoints += cells1D[index_1D_final[0]]

        #now the rest
        for k in range(1, len(index_1D_final)):
            aux = cells1D[index_1D_final[k]][1]
            indexToPoints.append(aux)
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #we will add now the scalars

        scalar_values = []  
                        #the value of the scalars pointed by elem in indexToPoints.
                        #scalar_values = [[temperature.values related to the nodes],
                        #                [dens.values related to the nodes], 
                        #                ...]

        for t in range(NumberOfTimeSlices):
            aux = []
            for scalar in OurScalars[t]:
                aux.append([scalar[node] for node in indexToPoints])
            scalar_values.append(aux)
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #Now we're going to take the coordinates of the nodes and make a list from
        #their distance, starting in 0 and ending in the correspondent distance.
        #this will be a great improvement for the 3D plot.
        
        def distance2D(p1, p2):
        #we will calculate distances with this function.    
            return  ((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)**0.5

        distances = [0] #this will be our list (firs element in position 0)

        for i in range(1, len(indexToPoints)):
            distances.append(
                distance2D(points_geo[i], points_geo[i-1])+distances[-1])
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #we write everything now



        outputName = "./temp_plugin/tunnel_plot.txt"
        g = open(outputName,"w")


        #first some information about the mesh
        g.write(chosen_name)            #the name of the subset
        g.write('\n')
        for position in distances:     #list for "x" variable for the plot
            g.write(str(position) + ' ')
        g.write('\n')

        for time in timeArray:          #list for "time" varibale
            g.write(str(time) + ' ')
        g.write('\n')


        #now the titles: "time", "point_ID" and scalars
        g.write("time" + ' ')
        g.write("point_ID" + ' ')
        for name in scalar_names:
            g.write(name + ' ')
        g.write("\n")


        print("writing everything on "+outputName+"...\n")


        for t in range(NumberOfTimeSlices):             #for every timestep
            for i in range(len(indexToPoints)):         #for every point
                g.write(str(t)+' ')
                g.write(str(indexToPoints[i])+' ')
                for elem in scalar_values[t]:           #for every scalar
                    g.write(str(elem[i]))
                    g.write(" ")
                g.write('\n')

        g.close()

        print("done!\n\n")
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    print("-----------------------------------------------------------------------")



    if TrackPoint:  #if we have decided to track a single point
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #we get here all the needed scalars

        print("preparing the data for PointTracker...")
        outputName = "./temp_plugin/trackpoint.txt"
        f = open(outputName,"w")


        OurScalars = [] #value = OurScalars[timestep][scalar (temp, dens...)]
        #this list contains the values of the sccalars for each timestep
        for t in range(NumberOfTimeSlices):
            scalar_names = [] 

            SCALARS, SC_INDEX, NAMES = ids_utils.ids_read_Scalars(t = t) #first t = 0
            aux = []
            for i, elem in enumerate(SCALARS):
                if len(elem) > 0:
                    scalar_names.append(NAMES[i])
                    aux.append(elem[ids_utils.GetScIndexForNodes()][point_ID])
            OurScalars.append(aux)
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #we write everything now

        f.write("time" + ' ')
        for i in range(len(scalar_names)):
            f.write(scalar_names[i] + ' ')
        f.write("\n")

        print("writing everything on "+outputName+"...\n")

        for t in range(NumberOfTimeSlices):
            
            f.write(str(timeArray[t]))
            f.write(" ")

            for elem in OurScalars[t]:
                f.write(str(elem))
                f.write(" ")
            f.write('\n')   
        #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        f.close()


    ids_utils.ids_close()
    print("-----------------------------------------------------------------------")
    print("Finished!")


#uncomment if running on a ParaView's Programmable Source

shot = 22108  
run = 1000
user = 'g2mcarpi'                       
device = 'imas20'           
version = '3'              
PrepareTrackSubset = True      #if you want to write a file for tracking a subset
PrepareTrackPoint = True       #if you want to track a point   
Subset = 5              #the subset that you want to track
Point = 900          #the point you want to track

CreateAndOpen(shot, run, user, device, version, Subset, \
                Point, PrepareTrackPoint, PrepareTrackSubset)

