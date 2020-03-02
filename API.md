# WIR3 API
**Document of all the APIs used in this project**

## Dataset
The class for each input dataset, it will keep track of the steps user want to perform and then execute them in order.
- **username**: name of the user
- **useremail**: email of the user
- **init()**: load the dataset into anndata through load function; document username and email
- **runner**: a runner object to control the precedure of the process and analysis
- **get_response()**: it will return the current processed anndata
- **analysis_steps()**: store the steps into a linked list with a state initialized to “unfinished” for each step, and a time intialized to 0 to record the processing time for each step; when user hit the submit button.

## Runner:
- **start()**: initialize this runner when a dataset is loaded
- **current_step()**: return the current step of the dataset
- **excute_function(api, state)**:it will call the function of current step api; set the state to “finish” when successfully execute without error.
- **get_response()**: it will call dataset’s get_response() to get the current data; and also it will return the current_step

## Proprocess:
The proprocess module with function we are trying to integrate:
- **method1()**:
- **method2()**:
- **method3()**:


## Analysis: The analysis module with function we are trying to integrate:
- **method1()**:
- **method2()**:
- **method3()**:

## Visualization: The analysis module corresponding to different visualization methods
- **method1()**:
- **method2()**:
- **method3()**:
