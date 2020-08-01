var jmol;
var Info = {
  color: "0xFFFFFF",
  height: 400,
  width: 550,
  use: "HTML5",
  allowJavaScript: false,
  j2sPath: "/static/jsmol/j2s",
  script: "set displaycellparameters false; load struct/cif { 1 1 1 };",
  addSelectionOptions: false,
  disableJ2SLoadMonitor: true,
  disableInitialConsole: true,
  debug: false,
};
var jmol2;
var Info2 = {
  color: "0xFFFFFF",
  height: 400,
  width: 550,
  use: "HTML5",
  allowJavaScript: false,
  j2sPath: "/static/jsmol/j2s",
  script: "set displaycellparameters false; load struct2/cif { 1 1 1 };",
  addSelectionOptions: false,
  disableJ2SLoadMonitor: true,
  disableInitialConsole: true,
  debug: false,
};

function repeatCell(n1, n2, n3, num) {
  var s = '{ ' + n1.toString() + ' ' + n2.toString() + ' ' + n3.toString() + ' };';
  if (num) Jmol.script(jmol2, 'load "" ' + s);
  else Jmol.script(jmol, 'load "" ' + s);
}

$(document).ready(function() {
  $("#appdiv").html(Jmol.getAppletHtml("jmol", Info));
  $("#appdiv2").html(Jmol.getAppletHtml("jmol2", Info2));
})
