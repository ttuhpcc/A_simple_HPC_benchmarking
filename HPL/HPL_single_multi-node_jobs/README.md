This contains files for running a single parallel job spanning **N** number of nodes to test the combined performance of the whole cluster.
1. First start with creating a **nodelist** containing a list of target nodes to run HPL.
    a. For cases with root access, use `create_node_list_snodes.sh` to create a nodelist. Change parameters accordingly.
    b. Else, use `create_node_list.sh`
2. Run `control_nodelist.sh` to read nodelist from previous step's output file.
    a. It also runs the `sample_script.sh` with a node argument which   
