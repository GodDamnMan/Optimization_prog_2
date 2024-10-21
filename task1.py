import math

spacing = 6

class Simplex:
    
    # we define tableu:list (n:int * m:int) - matrix for simplex
    # constrains:list - to check whether tableu ratios are in the scope of them
    # var - ['x1', 'x2', 'x3', 's1', 's2', 's3']
    # basic - ['s1', 's2', 's3']
    def __init__(self, obj_function: list, constrains_matrix: list, right_hand_side_num: list, epsilon:int):
        self.constrains = constrains_matrix
        self.tableu = [[-c for c in obj_function] + [0 for i in range(len(constrains_matrix))] + [0]]
        for i in range(len(constrains_matrix)):
            self.tableu.append([el for el in constrains_matrix[i]] + [1 if i == j else 0 for j in range(len(constrains_matrix))] + [right_hand_side_num[i]])
        
        self.n = len(self.tableu)
        self.m = len(self.tableu[0])
        self.eps = epsilon
        self.basic = [f's{i+1}' for i in range(len(constrains_matrix))]
        self.vars = [f'x{i+1}' for i in range(len(obj_function))] + self.basic
        self.formatting = '{:'+str(self.eps + spacing) + '.' + str(self.eps) + 'f}'
        self.solving = []

    # just simplex 
    def simplex_method(self):
        while True:
            enters = self.solving[0].index(min(self.solving[0]))
            if self.solving[0][enters] >= 0:
                break
            leaves = 0
            l_value = math.inf
            for i in range(1, self.n):
                if self.solving[i][enters] == 0: continue
                temp = self.solving[i][self.m-1]/self.solving[i][enters]
                if (temp < l_value):
                    leaves = i
                    l_value = temp
            if leaves == 0:
                break

            self.basic[leaves - 1] = self.vars[enters]
            
            for i in range(self.m):
                self.solving[leaves][i] /= self.solving[leaves][enters]
        
            for i in range(self.n):
                if i == leaves:
                    continue
                coef = -self.solving[i][enters] / self.solving[leaves][enters]
                for j in range(self.m):
                    self.solving[i][j] += self.solving[leaves][j] * coef

    
    def solve_maximize(self):
        self.solving = self.tableu.copy()
        self.simplex_method()
        return self.solving[0][self.m-1]
        
    
    def solve_minimize(self): 
        self.solving = self.tableu.copy()
        for i in range(self.m):
            self.solving[0][i] *= -1
        self.simplex_method()
        return -self.solving[0][self.m-1]

    
    # function to print initial tableu
    def print_initial(self):
        print("initial tableu:")
        string = '____'
        if ((self.eps + spacing)%2 == 1 or (self.eps + spacing) == 5):
            for j in range (len(self.tableu[0])-1):
                string += "|" + ((self.eps + spacing)//2) * " " + self.vars[j] + ((self.eps + spacing)//2) * " "
            string += "|" + ((self.eps + spacing)//2-1) * " " + "Sol" + (self.eps + spacing)//2 * " " + "|"
            print(string)
        else: 
            for j in range (len(self.tableu[0])-1):
                string += "|" + (self.eps + spacing)//2 * " " + self.vars[j] + ((self.eps + spacing)//2-1) * " "
            string += "|" + ((self.eps + spacing)//2-1) * " " + "Sol" + ((self.eps + spacing)//2-1) * " " + "|"
            print(string)
        k = 0
        for i in self.tableu:
            if(k == 0):
                print(" z  |", end = "")
            else:
                print(self.basic[k-1] + "  |", end = "")
            for j in i:
                print(self.formatting.format(j), end = " |")
            print()
            k+=1
        print()

    # function to print tableu after applying a Simplex method
    def print_solved(self):
        print("optimum is", self.solving[0][self.m-1])
        for i in range(len(self.basic)):
            print(self.basic[i], "=", self.solving[i+1][self.m-1])
            
        string = '____'
        if ((self.eps + 4)%2 == 1 or (self.eps + 4) == 5):
            for j in range (len(self.tableu[0])-1):
                string += "|" + ((self.eps + spacing)//2) * " " + self.vars[j] + ((self.eps + spacing)//2) * " "
            string += "|" + ((self.eps + spacing)//2-1) * " " + "Sol" + (self.eps + spacing)//2 * " " + "|"
            print(string)
        else: 
            for j in range (len(self.tableu[0])-1):
                string += "|" + (self.eps + spacing)//2 * " " + self.vars[j] + ((self.eps + spacing)//2-1) * " "
            string += "|" + ((self.eps + spacing)//2-1) * " " + "Sol" + ((self.eps + spacing)//2-1) * " " + "|"
            print(string)

        k = 0
        for i in self.solving:
            if(k == 0):
                print(" z  |", end = "")
            else:
                print(self.basic[k-1] + "  |", end = "")
            for j in i:
                print(self.formatting.format(j), end = " |")
            print()
            k+=1
        print()


def simplex_input():

    type = input("Greetings, this programm will solve your LP problem using Simplex method.\nEnter the type of the problem(Max/Min): ").lower()
    if (type != "max" and type != "min"):
        print("ERROR: UNKNOWN TYPE")
        return
    
    try:
        objective_function = list(map(float, input("Enter the coefficients of the objective function: ").split(" ")))
        if(len(objective_function) == 0):
            print("ERROR: NO COEFFICIENTS")
            return
        
        amount = int(input("Enter amount of the constraints(not assuming x>=0): "))
        if(amount < 1):
            print("ERROR: AMOUNT < 1 ?!")
        constraints = []
        print("Write constraints in format: num num sign\nExample: 3 4 <=")
        for i in range(amount):
            constraint = list(map(float, input(f"Enter the {i+1} constraint function coefficients: ").split(" ")))
            if(len(constraint) == 0):
                print("ERROR: NO COEFFICIENTS")
                return
            constraints.append(constraint)


        right_hand_side = list(map(float, input("Enter the right-hand side numbers: ").split(" ")))
        if(len(right_hand_side) != amount):
            print("ERROR: NOT ENOUGH COEFFICIENTS")
            return
        else:
            for i in right_hand_side:
                if(i < 0):
                    print("The method is not applicable!")
                    return

        accuracy = int(input("Enter the approximation accuracy: "))

        lp = Simplex(objective_function, constraints, right_hand_side, accuracy)
        lp.print_initial()
        if(type == "max"):
            lp.solve_maximize()
        else:
            lp.solve_minimize()
        lp.print_solved()


    except ValueError:
        print("ERROR: NOT A NUMBER")
        return
    