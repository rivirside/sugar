declare module "react-cytoscapejs" {
  import type cytoscape from "cytoscape";
  import type { ComponentType } from "react";

  interface CytoscapeComponentProps {
    id?: string;
    cy?: (cy: cytoscape.Core) => void;
    style?: React.CSSProperties;
    elements: cytoscape.ElementDefinition[];
    layout?: cytoscape.LayoutOptions;
    stylesheet?: cytoscape.StylesheetStyle[];
    className?: string;
    zoom?: number;
    pan?: cytoscape.Position;
    minZoom?: number;
    maxZoom?: number;
    zoomingEnabled?: boolean;
    userZoomingEnabled?: boolean;
    panningEnabled?: boolean;
    userPanningEnabled?: boolean;
    boxSelectionEnabled?: boolean;
    autoungrabify?: boolean;
    autounselectify?: boolean;
  }

  const CytoscapeComponent: ComponentType<CytoscapeComponentProps>;
  export default CytoscapeComponent;
}
