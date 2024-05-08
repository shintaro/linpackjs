var worker = new Array();

worker[0] = new Worker("linpackWorker.js");
worker[1] = new Worker("linpackWorker-compiled.js");
worker[2] = new Worker("linpackWorker-advanced-compiled.js");

worker[0].onerror = function(error) {
      alert(JSON.stringify(error));
};
worker[1].onerror = function(error) {
      alert(JSON.stringify(error));
};
worker[2].onerror = function(error) {
      alert(JSON.stringify(error));
};

window.onload = function(event) {
  document.getElementById("statusField").innerHTML = "Loaded";

  if (window.Worker) {
    document.getElementById("runButton").addEventListener("click", send, false);
  } else {
    alert("This browser doesn't support Web Wokers");
  }  
}

function send() { 
  var arraySize = document.getElementById("arraySize").value;
  document.getElementById("resultField").innerHTML = "";
  if (arraySize > 99 && arraySize < 5001) {
      var code = document.getElementsByName("code");
      for(var i = 0; i < code.length; i++){
          if (code[i].checked) {
              worker[i].postMessage({ "arraySize": arraySize });
          }
      }
      document.getElementById("statusField").innerHTML = "Running";
  } else {
      alert("Please enter valid number");
  }
}

worker[0].onmessage = function (event) {
  document.getElementById("statusField").innerHTML = "Ready";
  printResult(event.data);
}
worker[1].onmessage = function (event) {
  document.getElementById("statusField").innerHTML = "Ready";
  printResult(event.data);
}
worker[2].onmessage = function (event) {
  document.getElementById("statusField").innerHTML = "Ready";
  printResult(event.data);
}

function printResult(s) {
  var str = s.split(" ");
  var text = "LINPACK benchmark, Double precision.<br>";
  text += "Array size ";
  text += document.getElementById("arraySize").value; 
  text += " X ";
  text += document.getElementById("arraySize").value; 
  text += "<br><br>";
  text += "Average rolled and unrolled performance:";
  text += "<table border='1'>";
  text += "<tr>";
  text += "<td>Reps</td><td>Time(s)</td><td>DGEFA</td><td>DGESL</td><td>OVERHEAD</td><td>KFOPS</td>";
  text += "</tr>";
  text += "<tr>";
  text += "<td>" + str[0] + "</td>";
  text += "<td>" + str[1] + "</td>";
  text += "<td>" + str[2] + "</td>";
  text += "<td>" + str[3] + "</td>";
  text += "<td>" + str[4] + "</td>";
  text += "<td>" + str[5] + "</td>";
  text += "</tr>";
  text += "</table>";
  document.getElementById("resultField").innerHTML = text;
}

